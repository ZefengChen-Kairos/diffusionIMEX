"""Reference silhouettes from the videos and a simulation-to-reference similarity score.

The cup surface in a video frame is an ellipse (oblique view).  We segment the liquid (crema +
foam inside the cup), fit its ellipse from the second moments of the filled convex hull, map it
affinely onto the unit disc, and segment the white foam there.  The similarity of a simulated
pattern is the white-mask IoU on the unit disc, maximised over rotations and a mirror (the camera
orientation of each video is arbitrary).

    python -m gates.reference            # writes reference masks + previews to <out>/reference/
"""
import os
import json
import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.spatial import ConvexHull
from matplotlib.path import Path

REF = {
    # name: (video file suffix, time [s])
    "heart": ("reference_heart.mp4", 11.3),
    "push_heart": ("reference_push_heart.mp4", 15.5),
    "layered_heart": ("reference_layered_heart.mp4", 11.0),
    "tulip": ("reference_layered_heart.mp4", 12.5),     # top view "No.4 压纹郁金香"
    "leaf": ("reference_leaf.mp4", 13.5),
    "swan": ("reference_leaf.mp4", 15.3),               # top view "No.6 天鹅"
}
NU = 256          # unit-disc raster


def _hsv(a):
    a = a.astype(float) / 255.0
    mx = a.max(-1)
    mn = a.min(-1)
    sat = (mx - mn) / np.maximum(mx, 1e-6)
    return mx, sat


def liquid_mask(img):
    """Crema (saturated brown) region of the cup, filled by its convex hull."""
    v, s = _hsv(img)
    r, g, b = [img[..., k].astype(float) for k in range(3)]
    crema = (s > 0.30) & (v > 0.35) & (r > g) & (g > b) & (r - b > 50)
    lab, n = ndimage.label(crema)
    if n == 0:
        raise ValueError("no crema found")
    sizes = ndimage.sum(crema, lab, range(1, n + 1))
    big = lab == (1 + int(np.argmax(sizes)))
    ys, xs = np.nonzero(big)
    pts = np.stack([xs, ys], 1).astype(float)
    hull = ConvexHull(pts)
    path = Path(pts[hull.vertices])
    H, W = big.shape
    yy, xx = np.mgrid[0:H, 0:W]
    inside = path.contains_points(np.stack([xx.ravel(), yy.ravel()], 1)).reshape(H, W)
    return inside


def ellipse_map(mask):
    """Centre and affine map image-pixel -> unit disc from the filled-ellipse moments."""
    ys, xs = np.nonzero(mask)
    c = np.array([xs.mean(), ys.mean()])
    X = np.stack([xs - c[0], ys - c[1]], 1)
    C = X.T @ X / len(X)
    w, V = np.linalg.eigh(C)
    A = V @ np.diag(1.0 / (2.0 * np.sqrt(w))) @ V.T      # unit disc: A (x - c)
    return c, A


def to_unit_disc(img, c, A, n=NU, shrink=0.93):
    """Resample the image on the unit disc (u, v in [-1, 1]); v points up in the image."""
    Ainv = np.linalg.inv(A)
    u = np.linspace(-1, 1, n)
    U, Vv = np.meshgrid(u, -u)
    pts = np.stack([U.ravel(), Vv.ravel()], 0) * shrink
    xy = (Ainv @ pts) + c[:, None]
    out = np.zeros((n * n, 3))
    for k in range(3):
        out[:, k] = ndimage.map_coordinates(img[..., k].astype(float), [xy[1], xy[0]], order=1, mode="nearest")
    return out.reshape(n, n, 3), (U ** 2 + Vv ** 2) < 1.0


def _otsu(x):
    h, e = np.histogram(x, bins=128)
    c = 0.5 * (e[1:] + e[:-1])
    w0 = np.cumsum(h)
    w1 = w0[-1] - w0
    m0 = np.cumsum(h * c) / np.maximum(w0, 1)
    m1 = ((h * c).sum() - np.cumsum(h * c)) / np.maximum(w1, 1)
    return c[int(np.argmax(w0 * w1 * (m0 - m1) ** 2))]


def white_mask_ref(disc_img, inside):
    """Foam = bright and unsaturated; Otsu threshold of the whiteness v (1 - s) inside the cup."""
    v, s = _hsv(disc_img)
    w = v * (1.0 - s)
    thr = _otsu(w[inside])
    return (w > thr) & inside


def white_mask_sim(alpha, grid, n=NU, shrink=0.93):
    """Resample a simulated alpha (fine grid, y up = forward) onto the unit disc (same shrink as the
    reference: the unit disc is shrink x the liquid radius)."""
    u = np.linspace(-1, 1, n)
    U, Vv = np.meshgrid(u, -u)              # image convention: row 0 is the top
    x = U * grid.R * shrink
    y = Vv * grid.R * shrink
    j = (y + 0.5 * grid.L) / grid.dx - 0.5
    i = (x + 0.5 * grid.L) / grid.dx - 0.5
    a = ndimage.map_coordinates(alpha, [j, i], order=1, mode="nearest")
    return (a > 0.5) & ((U ** 2 + Vv ** 2) < 1.0)


def best_iou(sim_mask, ref_mask, step=5):
    """Max IoU over rotations (step degrees) and a mirror; returns (iou, angle, mirrored, aligned mask)."""
    best = (-1.0, 0, False, None)
    for mir in (False, True):
        m = sim_mask[:, ::-1] if mir else sim_mask
        for ang in range(0, 360, step):
            r = ndimage.rotate(m.astype(float), ang, reshape=False, order=1) > 0.5
            inter = (r & ref_mask).sum()
            uni = (r | ref_mask).sum()
            iou = inter / uni if uni else 0.0
            if iou > best[0]:
                best = (float(iou), ang, mir, r)
    return best


def video_path(suffix):
    up = "/root/.claude/uploads"
    for root, _, files in os.walk(up):
        for f in files:
            if f.endswith(suffix):
                return os.path.join(root, f)
    local = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reference", suffix)
    return local if os.path.exists(local) else None


def frame(path, t):
    import subprocess
    import tempfile
    tmp = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
    subprocess.run(["ffmpeg", "-v", "error", "-y", "-ss", str(t), "-i", path, "-frames:v", "1", tmp], check=True)
    return np.asarray(Image.open(tmp).convert("RGB"))


def build(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    res = {}
    for name, (suffix, t) in REF.items():
        p = video_path(suffix)
        png = os.path.join(out_dir, f"{name}_frame.png")
        if p is not None:
            img = frame(p, t)
            Image.fromarray(img).save(png)
        elif os.path.exists(png):
            img = np.asarray(Image.open(png).convert("RGB"))
        else:
            continue
        liq = liquid_mask(img)
        c, A = ellipse_map(liq)
        disc, inside = to_unit_disc(img, c, A)
        wm = white_mask_ref(disc, inside)
        wm = ndimage.binary_opening(wm, iterations=1)
        np.save(os.path.join(out_dir, f"{name}_mask.npy"), wm)
        prev = np.concatenate([disc.astype(np.uint8), np.repeat((wm * 255).astype(np.uint8)[..., None], 3, 2)], 1)
        Image.fromarray(prev).save(os.path.join(out_dir, f"{name}_disc.png"))
        res[name] = dict(white_fraction=float(wm.sum() / inside.sum()))
    json.dump(res, open(os.path.join(out_dir, "reference.json"), "w"), indent=1)
    return res


def load_masks(out_dir):
    m = {}
    for name in REF:
        p = os.path.join(out_dir, f"{name}_mask.npy")
        if os.path.exists(p):
            m[name] = np.load(p)
    return m


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from gates.common import out_root
    print(build(os.path.join(out_root(), "reference")))


def score_runs(runs, out_dir, P=None, NP=None, tag="similarity"):
    """runs: {name: (Yf, M, sim_or_shell)} -> IoU table + aligned overlay board."""
    from gates.common import alpha_of
    from render.render import board
    from gates.common import out_root
    refs = load_masks(os.path.join(out_root(), "reference"))
    res, rows = {}, []
    for name, (Y, M, sh) in runs.items():
        if name not in refs:
            continue
        a = alpha_of(sh, Y, M)
        sm = white_mask_sim(a, sh.gf)
        iou, ang, mir, al = best_iou(sm, refs[name])
        res[name] = dict(iou=iou, angle=ang, mirror=mir, white_sim=float(sm.sum() / (np.pi * (NU / 2) ** 2)),
                         white_ref=float(refs[name].sum() / (np.pi * (NU / 2) ** 2)))
        ov = np.zeros(al.shape + (3,), np.uint8) + 30
        ov[refs[name]] = (200, 120, 60)
        ov[al] = (90, 160, 230)
        ov[al & refs[name]] = (245, 245, 245)
        rows.append(ov)
    if rows:
        board([rows], os.path.join(out_dir, f"{tag}.png"), labels=[f"{n} IoU {res[n]['iou']:.2f}" for n in res],
              cell=220)
    return res
