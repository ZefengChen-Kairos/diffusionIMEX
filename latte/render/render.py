"""Optical rendering, separated from the dynamics (section 1.3 E).

alpha = 1 - exp(-(Phi/Phi_*)^2);  I = (1 - alpha) I_crema + alpha I_milk.
One shared Phi_*, one crema colour with a fixed radial gradient, one milk colour.
"""
import numpy as np
from PIL import Image, ImageDraw, ImageFont

# colours sampled from the reference video frames (median crema / rim crema / foam)
I_CREMA_C = np.array([208, 158, 100], float)   # crema near the centre
I_CREMA_E = np.array([150, 92, 56], float)     # crema at the rim (radial gradient)
I_RIM = np.array([112, 64, 38], float)         # meniscus line at the wall
I_MILK = np.array([251, 248, 241], float)
I_CUP = np.array([243, 242, 238], float)
I_BG = np.array([255, 255, 255], float)
SS = 2                                          # supersampling of the fine grid for display


def alpha(phi, Phi_star):
    return 1.0 - np.exp(-(phi / Phi_star) ** 2)


def _smooth(x, a, b):
    t = np.clip((x - a) / (b - a), 0.0, 1.0)
    return t * t * (3 - 2 * t)


def render(phi, grid, Phi_star, flip=True, ss=SS):
    """phi on the fine grid -> uint8 RGB image (y up when flip=True).

    Same formula for every pattern: alpha = 1 - exp(-(Phi/Phi_*)^2), I = (1-alpha) I_crema + alpha I_milk,
    with a fixed radial crema gradient.  For display, Phi is upsampled (cubic) by `ss` before alpha
    and the cup wall / rim are drawn around the liquid."""
    from scipy import ndimage
    if ss > 1:
        phi = np.maximum(ndimage.zoom(phi, ss, order=3, mode="nearest"), 0.0)
    n = phi.shape[0]
    xc = (np.arange(n) - 0.5 * (n - 1)) * (grid.L / n)
    X, Y = np.meshgrid(xc, xc)
    rr = np.sqrt(X ** 2 + Y ** 2) / grid.R
    a = alpha(phi, Phi_star)
    a = a * (rr < 1.0)
    a = a[..., None]
    w = _smooth(rr, 0.35, 1.0)[..., None] ** 1.5
    crema = (1 - w) * I_CREMA_C + w * I_CREMA_E
    m = _smooth(rr, 0.965, 0.995)[..., None]
    crema = (1 - m) * crema + m * I_RIM
    img = (1 - a) * crema + a * I_MILK
    # cup wall: liquid edge anti-aliased over ~1 px, ceramic rim, white background
    px = grid.L / n / grid.R
    cov = np.clip((1.0 - rr) / px + 0.5, 0.0, 1.0)[..., None]
    rim = _smooth(rr, 1.0, 1.02)[..., None] * (1 - _smooth(rr, 1.10, 1.115)[..., None])
    shade = (0.93 + 0.07 * _smooth(rr, 1.0, 1.1))[..., None]
    cup = I_CUP * shade
    outer = rim * cup + (1 - rim) * I_BG
    outer = np.where((rr < 1.05)[..., None], cup, outer)
    img = cov * img + (1 - cov) * outer
    img = np.clip(img, 0, 255).astype(np.uint8)
    return img[::-1] if flip else img


def tag_image(Yf, grid, flip=True):
    """False-colour image of the batch tags (diagnostic only)."""
    pal = np.array([[230, 60, 60], [60, 160, 230], [80, 200, 90], [240, 190, 40],
                    [170, 90, 220], [240, 120, 200], [120, 120, 120]], float)
    K = Yf.shape[0] - 1
    tot = np.maximum(Yf[0], 1e-12)
    col = np.zeros(Yf.shape[1:] + (3,))
    for k in range(K):
        col += (Yf[1 + k] / tot)[..., None] * pal[k % len(pal)]
    a = np.clip(Yf[0] / 0.5, 0, 1)[..., None]
    img = (1 - a) * np.array([60, 40, 30.0]) + a * col
    img = np.where(grid.wet[..., None], img, I_CUP)
    img = np.clip(img, 0, 255).astype(np.uint8)
    return img[::-1] if flip else img


def save(img, path, scale=1):
    im = Image.fromarray(img)
    if scale != 1:
        im = im.resize((im.width * scale, im.height * scale), Image.LANCZOS)
    im.save(path)


def board(rows, path, labels=None, row_labels=None, cell=220, pad=6):
    """rows: list of lists of uint8 images (or None)."""
    nr = len(rows)
    nc = max(len(r) for r in rows)
    top = 22 if labels else 0
    left = 110 if row_labels else 0
    W = left + nc * (cell + pad) + pad
    H = top + nr * (cell + pad) + pad
    canvas = Image.new("RGB", (W, H), (255, 255, 255))
    d = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 13)
    except OSError:
        font = ImageFont.load_default()
    for r, row in enumerate(rows):
        for c, im in enumerate(row):
            if im is None:
                continue
            pim = Image.fromarray(im) if isinstance(im, np.ndarray) else im
            pim = pim.resize((cell, cell), Image.LANCZOS) if pim.size != (cell, cell) else pim
            canvas.paste(pim, (left + pad + c * (cell + pad), top + pad + r * (cell + pad)))
        if row_labels:
            d.text((6, top + pad + r * (cell + pad) + cell // 2 - 8), row_labels[r], fill=(0, 0, 0), font=font)
    if labels:
        for c, lb in enumerate(labels):
            d.text((left + pad + c * (cell + pad) + 4, 4), lb, fill=(0, 0, 0), font=font)
    canvas.save(path)
    return canvas


def crop_cup(img, grid, margin=0.0055):
    """Crop a rendered image to the cup disc bounding box."""
    N = img.shape[0]
    half = (grid.R + margin) / grid.L * N
    c = N / 2
    a, b = int(max(c - half, 0)), int(min(c + half, N))
    return img[a:b, a:b]
