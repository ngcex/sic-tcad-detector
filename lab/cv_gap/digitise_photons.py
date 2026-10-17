"""Digitise Photons Fig. 4b (eps*A/C, um, vs bias) from the native raster of the PDF.

usage: python digitise_photons.py fig4-000.png out.csv overlay.png
Panel (b) = right half. Axis calibration from the OUTWARD major tick marks on the
frame; the orange hollow-diamond markers (one per 250 mV step) carry the data.
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
# ---- frame of panel (b)
col_dark = dark[:, W // 2 :].sum(axis=0)
xs = np.where(col_dark > 0.6 * H)[0] + W // 2
row_dark = dark[:, W // 2 :].sum(axis=1)
ys = np.where(row_dark > 0.3 * W)[0]
x_left, x_right = xs.min(), xs.max()
y_top, y_bot = ys.min(), ys.max()
print(f"frame x {x_left}..{x_right}  y {y_top}..{y_bot}")


def runs(mask1d, offset):
    idx = np.where(mask1d)[0]
    if len(idx) == 0:
        return np.array([])
    groups = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
    return np.array([gr.mean() + offset for gr in groups])


# major x ticks: rows y_bot+7 .. y_bot+10 (outward, longer than minor)
xt = runs(
    dark[y_bot + 7 : y_bot + 11, x_left - 2 : x_right + 3].all(axis=0), x_left - 2
)
yt = runs(dark[y_top - 2 : y_bot + 3, x_left - 11 : x_left - 7].all(axis=1), y_top - 2)
print("x major ticks px:", np.round(xt, 1))
print("y major ticks px:", np.round(yt, 1))
V_ticks = np.arange(-60, 1, 10)[: len(xt)]
assert len(xt) == 7, len(xt)
# y: 0,2,...,10 from the bottom (the 0 tick may merge with the frame; detect count)
yt = np.sort(yt)[::-1]  # bottom (largest px) first
W_ticks = 2.0 * np.arange(len(yt))
ax_, bx_ = np.polyfit(xt, V_ticks, 1)
ay_, by_ = np.polyfit(yt, W_ticks, 1)
print(
    f"x: {1/ax_:.3f} px/V, resid {np.abs(V_ticks-(ax_*xt+bx_)).max():.3f} V; "
    f"y: {1/abs(ay_):.2f} px/um, resid {np.abs(W_ticks-(ay_*yt+by_)).max():.4f} um"
)
V_of = lambda px: ax_ * px + bx_
W_of = lambda py: ay_ * py + by_
px_of_V = lambda V: (V - bx_) / ax_

# ---- orange marker pixels (MATLAB orange 217,83,25; anti-aliased edges included)
orange = (r > 150) & (b < 140) & (r - b > 80) & (g < 170)
orange[:, : x_left + 2] = False
orange[y_bot - 1 :, :] = False
# weight: how orange (0..1)
wgt = np.clip((r - b) / 192.0, 0, 1) * orange

# --- A. dense band: per-column weighted centroid of orange, in a window around the
# local band (robust to the blue line, which is not orange)
V_grid = np.arange(-60.0, 2.5001, 0.25)
rows = np.arange(H)
res = []
for V in V_grid:
    xc = px_of_V(V)
    c0, c1 = int(np.floor(xc)), int(np.ceil(xc))
    fr = xc - c0
    col = (1 - fr) * wgt[:, c0] + fr * wgt[:, c1]
    if col.sum() < 1.0:
        res.append((V, np.nan, np.nan, 0.0))
        continue
    yc = (rows * col).sum() / col.sum()
    sd = np.sqrt(((rows - yc) ** 2 * col).sum() / col.sum())
    res.append((V, W_of(yc), sd * abs(ay_), col.sum()))
res = np.array(res)

# --- B. isolated markers (knee): connected components of the orange mask
lab, n = ndimage.label(orange)
sizes = ndimage.sum(orange, lab, range(1, n + 1))
cents = ndimage.center_of_mass(wgt, lab, range(1, n + 1))
markers = []
objs = ndimage.find_objects(lab)
for i, (sz, (cy, cx)) in enumerate(zip(sizes, cents)):
    sl = objs[i]
    hgt = sl[0].stop - sl[0].start
    wid = sl[1].stop - sl[1].start
    markers.append((cx, cy, sz, hgt, wid))
markers = np.array(markers)
print(f"{n} orange components; sizes median {np.median(sizes):.0f}")
single = markers[(markers[:, 3] < 20) & (markers[:, 4] < 20)]
print(f"{len(single)} isolated diamonds (h,w < 20 px):")
for cx, cy, sz, hgt, wid in sorted(single, key=lambda m: m[0]):
    print(
        f"   V={V_of(cx):7.3f}  W={W_of(cy):6.3f}  size {sz:.0f} h {hgt:.0f} w {wid:.0f}"
    )

# ---- write
with open(sys.argv[2], "w") as f:
    f.write("V,W_um,band_sd_um,weight\n")
    for V, Wv, sd, wt in res:
        f.write(f"{V:.2f},{Wv:.4f},{sd:.4f},{wt:.1f}\n")
with open(sys.argv[2].replace(".csv", "_markers.csv"), "w") as f:
    f.write("V,W_um,size,h,w\n")
    for cx, cy, sz, hgt, wid in sorted(single, key=lambda m: m[0]):
        f.write(f"{V_of(cx):.3f},{W_of(cy):.4f},{sz:.0f},{hgt:.0f},{wid:.0f}\n")

# ---- overlay
ov = Image.open(sys.argv[1]).convert("RGB")
d = ImageDraw.Draw(ov)
for V, Wv, sd, wt in res:
    if np.isfinite(Wv):
        x = px_of_V(V)
        y = (Wv - by_) / ay_
        d.ellipse([x - 2, y - 2, x + 2, y + 2], outline=(0, 160, 0))
for cx, cy, *_ in single:
    d.rectangle([cx - 3, cy - 3, cx + 3, cy + 3], outline=(160, 0, 160))
for t in xt:
    d.line([t, y_bot, t, y_bot + 30], fill=(0, 200, 0))
for t in yt:
    d.line([x_left - 30, t, x_left, t], fill=(0, 200, 0))
ov.crop((x_left - 60, y_top - 10, x_right + 10, y_bot + 40)).save(sys.argv[3])
