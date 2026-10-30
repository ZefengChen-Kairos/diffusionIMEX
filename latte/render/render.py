"""Optical rendering, separated from the dynamics (section 1.3 E).

alpha = 1 - exp(-(Phi/Phi_*)^2);  I = (1 - alpha) I_crema + alpha I_milk.
One shared Phi_*, one crema colour with a fixed radial gradient, one milk colour.
"""
import numpy as np
from PIL import Image, ImageDraw, ImageFont

I_CREMA_C = np.array([176, 112, 58], float)    # centre crema
I_CREMA_E = np.array([128, 72, 34], float)     # rim crema (radial gradient)
I_MILK = np.array([250, 244, 230], float)
I_CUP = np.array([236, 236, 232], float)


def alpha(phi, Phi_star):
    return 1.0 - np.exp(-(phi / Phi_star) ** 2)


def render(phi, grid, Phi_star, flip=True):
    """phi on the fine grid -> uint8 RGB image (y up when flip=True)."""
    a = alpha(phi, Phi_star)[..., None]
    rr = np.sqrt(grid.X ** 2 + grid.Y ** 2) / grid.R
    w = np.clip(rr, 0, 1)[..., None] ** 2
    crema = (1 - w) * I_CREMA_C + w * I_CREMA_E
    img = (1 - a) * crema + a * I_MILK
    img = np.where(grid.wet[..., None], img, I_CUP)
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


def crop_cup(img, grid, margin=0.004):
    """Crop a rendered image to the cup disc bounding box."""
    N = img.shape[0]
    half = (grid.R + margin) / grid.L * N
    c = N / 2
    a, b = int(max(c - half, 0)), int(min(c + half, N))
    return img[a:b, a:b]
