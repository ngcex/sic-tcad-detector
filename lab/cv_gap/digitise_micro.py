"""Digitise Microdosimetry Fig. 3b (Petringa 2025 Physica Medica 140, p.4): red filled
circles = 'Depletion layer [um]' (right axis 0..10), blue = capacitance (left axis,
labelled nF/cm^2 but the 1e-10 scale and the 4 mm2 plateau say F), 250 mV steps.
usage: python digitise_micro.py img-001.png out_prefix
Axes: bottom major ticks (-60..0 V, outward), left ticks (0.5..4.0 e-10), right red
tick LABELS (2..10 um; matplotlib centres the label on the tick).
"""

import sys

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage

img = np.array(Image.open(sys.argv[1]).convert("RGB")).astype(int)
H, W = img.shape[:2]
r, g, b = img[..., 0], img[..., 1], img[..., 2]
s = r + g + b
dark = s < 300
nonwhite = s < 600
xs = [x for x in range(W) if dark[:, x].sum() > 0.5 * H]
ys = [y for y in range(H) if dark[y, :].sum() > 0.5 * W]
x_left, x_right, y_top, y_bot = min(xs), max(xs), min(ys), max(ys)


def runs(mask1d, offset):
    idx = np.where(mask1d)[0]
    groups = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
    return np.array([gr.mean() + offset for gr in groups if len(gr)])


xt = runs(dark[y_bot + 3 : y_bot + 5, x_left - 2 : x_right + 3].all(axis=0), x_left - 2)
assert len(xt) == 7, xt
ax_, bx_ = np.polyfit(xt, np.arange(-60, 1, 10), 1)
yt_b = runs(
    nonwhite[y_top - 2 : y_bot + 3, x_left - 5 : x_left - 2].all(axis=1), y_top - 2
)
assert len(yt_b) == 8, yt_b
aC, bC = np.polyfit(yt_b, np.arange(4.0, 0.49, -0.5) * 1e-10, 1)
# right axis: red label centroids (rows), from the pixel dump (labels 10, 8, 6, 4, 2)
yt_r = np.array([35.0, 139.5, 244.5, 350.5, 455.5])
aW, bW = np.polyfit(yt_r, np.array([10.0, 8.0, 6.0, 4.0, 2.0]), 1)
print(
    f"x: {1 / ax_:.3f} px/V resid {np.abs(np.arange(-60, 1, 10) - (ax_ * xt + bx_)).max():.3f} V; "
    f"C: {abs(1 / aC) * 1e-10:.1f} px per 1e-10 F, resid {np.abs(np.arange(4.0, 0.49, -0.5) * 1e-10 - (aC * yt_b + bC)).max():.2e} F; "
    f"W: {abs(1 / aW):.2f} px/um, label resid {np.abs(np.array([10, 8, 6, 4, 2]) - (aW * yt_r + bW)).max():.3f} um"
)
V_of = lambda px: ax_ * px + bx_
px_of_V = lambda V: (V - bx_) / ax_
W_of = lambda py: aW * py + bW
C_of = lambda py: aC * py + bC

red = (r > 150) & (g < 110) & (b < 110)
blue = (b > 140) & (r < 110) & (g < 170)
for m in (red, blue):
    m[:, : x_left + 1] = False
    m[:, x_right:] = False
    m[y_bot:, :] = False
    m[: y_top + 1, :] = False
rows = np.arange(H)
V_grid = np.arange(-60.0, 2.5001, 0.25)


def band(mask, conv):
    out = []
    for V in V_grid:
        xc = px_of_V(V)
        c0 = int(np.floor(xc))
        fr = xc - c0
        col = (1 - fr) * mask[:, c0] + fr * mask[:, c0 + 1]
        if col.sum() < 1:
            out.append((V, np.nan, np.nan, 0))
            continue
        yc = (rows * col).sum() / col.sum()
        sd = np.sqrt(((rows - yc) ** 2 * col).sum() / col.sum())
        out.append((V, conv(yc), sd, col.sum()))
    return np.array(out)


Wb = band(red.astype(float), W_of)
Cb = band(blue.astype(float), C_of)
lab, n = ndimage.label(red)
objs = ndimage.find_objects(lab)
cents = ndimage.center_of_mass(red, lab, range(1, n + 1))
single = []
for i, (cy, cx) in enumerate(cents):
    sl = objs[i]
    h, w = sl[0].stop - sl[0].start, sl[1].stop - sl[1].start
    if 5 <= h <= 12 and w <= 12:
        single.append((V_of(cx), W_of(cy), h, w))
single = np.array(sorted(single))
print(f"{len(single)} isolated red markers")
for V, Wv, h, w in single:
    print(f"   V={V:7.3f}  W={Wv:6.3f}  ({h:.0f}x{w:.0f} px)")
with open(sys.argv[2] + "_band.csv", "w") as f:
    f.write("V,W_um,W_band_sd_px,W_weight,C_F,C_band_sd_px,C_weight\n")
    for (V, Wv, sdw, ww), (_, Cv, sdc, wc) in zip(Wb, Cb):
        f.write(f"{V:.2f},{Wv:.4f},{sdw:.2f},{ww:.0f},{Cv:.4e},{sdc:.2f},{wc:.0f}\n")
with open(sys.argv[2] + "_markers.csv", "w") as f:
    f.write("V,W_um,h,w\n")
    for V, Wv, h, w in single:
        f.write(f"{V:.3f},{Wv:.4f},{h:.0f},{w:.0f}\n")
ov = Image.open(sys.argv[1]).convert("RGB")
d = ImageDraw.Draw(ov)
for V, Wv, *_ in Wb:
    if np.isfinite(Wv):
        x, y = px_of_V(V), (Wv - bW) / aW
        d.ellipse([x - 1, y - 1, x + 1, y + 1], outline=(0, 160, 0))
for V, Cv, *_ in Cb:
    if np.isfinite(Cv):
        x, y = px_of_V(V), (Cv - bC) / aC
        d.ellipse([x - 1, y - 1, x + 1, y + 1], outline=(0, 160, 0))
for V, Wv, *_ in single:
    x, y = px_of_V(V), (Wv - bW) / aW
    d.rectangle([x - 3, y - 3, x + 3, y + 3], outline=(160, 0, 160))
ov.save(sys.argv[2] + "_overlay.png")
