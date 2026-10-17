"""Digitise Petringa 2025 (Physica Medica 140) Fig. 4b, blue CCE markers, from the PDF.

usage: uv run --with pillow python lab/microdosimeter/bench/digitise_fig4b.py "<path>/ Microdosimetry.pdf"

Renders page 4 at 400 dpi (pdftoppm), crops Fig. 4b, isolates the blue markers by
eroding away the 4-px connecting line (disk r=6), and maps marker centroids through the
outer y-axis major ticks (1.00 ... 0.92) and the x positions of the -50 and 0 V markers.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

pdf = sys.argv[1]
with tempfile.TemporaryDirectory() as d:
    subprocess.run(
        ["pdftoppm", "-f", "4", "-l", "4", "-r", "400", "-png", pdf, f"{d}/p"],
        check=True,
    )
    img = Image.open(next(Path(d).glob("p*.png"))).convert("RGB")
w, h = img.size
im = np.array(
    img.crop((int(w * 0.49), int(h * 0.30), int(w * 0.95), int(h * 0.55)))
).astype(int)
r, g, b = im[..., 0], im[..., 1], im[..., 2]
blue = (b > 150) & (r < 80) & (g > 90) & (g < 150)
yy, xx = np.mgrid[-6:7, -6:7]
lab, n = ndimage.label(ndimage.binary_erosion(blue, (xx**2 + yy**2) <= 36))
markers = sorted(
    (np.where(lab == i)[1].mean(), np.where(lab == i)[0].mean())
    for i in range(1, n + 1)
)
# outer major ticks of the left (CCE) axis: dark pixels in column ~153 of the crop
nonwhite = (r + g + b) < 600
col = nonwhite[:, 153]
rows = np.where(col)[0]
groups = np.split(rows, np.where(np.diff(rows) > 2)[0] + 1)
ticks = np.array([gr.mean() for gr in groups if 200 < gr.mean() < 900])
labels = np.array([1.00, 0.98, 0.96, 0.94, 0.92])[: len(ticks)]
a, c = np.polyfit(ticks, labels, 1)
resid = labels - (a * ticks + c)
print(
    f"ticks px {np.round(ticks, 1)}  ({1 / abs(a):.0f} px per unit CCE), fit resid {np.abs(resid).max():.1e}"
)
xm = np.array([m[0] for m in markers])
volts = -50.0 + (xm - xm[0]) / (xm[-1] - xm[0]) * 50.0
for (xp, yp), V in zip(markers, volts):
    print(f"V = {V:6.2f}   CCE = {a * yp + c:.4f}   (marker px {xp:.1f}, {yp:.1f})")
