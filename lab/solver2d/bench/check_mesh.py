"""STEP 0 checks on a dumped devsim mesh npz (tensor grid, diagonal couplings, volumes)."""
import sys
import numpy as np

d = np.load(sys.argv[1])
x, y, vol, ec, E = d["x"], d["y"], d["NodeVolume"], d["EdgeCouple"], d["edges"]
xs, ys = np.unique(x), np.unique(y)
print("nodes", len(x), "nx", len(xs), "ny", len(ys), "nx*ny", len(xs) * len(ys))
print("x grid (um)", np.round(xs * 1e4, 4))
i = np.searchsorted(xs, x); j = np.searchsorted(ys, y)
assert np.all(xs[i] == x) and np.all(ys[j] == y)
key = j * len(xs) + i
print("unique (i,j):", len(np.unique(key)) == len(x))
di = np.abs(i[E[:, 0]] - i[E[:, 1]]); dj = np.abs(j[E[:, 0]] - j[E[:, 1]])
hor = (di == 1) & (dj == 0); ver = (di == 0) & (dj == 1); dia = (di == 1) & (dj == 1)
print("edges hor/ver/diag/other", hor.sum(), ver.sum(), dia.sum(), (~(hor | ver | dia)).sum())
print("diag |EdgeCouple| max", np.abs(ec[dia]).max(), " min hor/ver couple", ec[hor].min(), ec[ver].min())
# orientation of diagonals: does (i,j)->(i+1,j+1) or (i+1,j)->(i,j+1)?
s = (i[E[dia, 1]] - i[E[dia, 0]]) * (j[E[dia, 1]] - j[E[dia, 0]])
print("diag orientation sign counts (+1 = /, -1 = \\)", np.unique(s, return_counts=True))
# rectangular dual: ec(horizontal) == dy_dual, ec(vertical) == dx_dual
hx = np.diff(xs); hy = np.diff(ys)
dxd = np.zeros(len(xs)); dxd[:-1] += hx / 2; dxd[1:] += hx / 2
dyd = np.zeros(len(ys)); dyd[:-1] += hy / 2; dyd[1:] += hy / 2
vol_rect = dxd[i] * dyd[j]
print("NodeVolume vs rect dual: max rel", np.max(np.abs(vol - vol_rect) / vol_rect))
jh = j[E[hor, 0]]; ih = i[E[ver, 0]]
print("EdgeCouple hor vs dy_dual max rel", np.max(np.abs(ec[hor] - dyd[jh]) / dyd[jh]))
print("EdgeCouple ver vs dx_dual max rel", np.max(np.abs(ec[ver] - dxd[ih]) / dxd[ih]))
print("air buffer: y min/max", ys[0], ys[-1], "contact rows y=0 idx", np.where(ys == 0)[0], " y=max idx", len(ys) - 1)
print("contact lists", d["contact_anode"][:5], d["contact_cathode"][:5])
print("node order: first 5 (i,j)", list(zip(i[:5], j[:5])), " is j fastest?", np.all(key[np.argsort(key)] == key[np.argsort(key)]))
print("devsim index -> j*nx+i identity?", np.array_equal(key, np.arange(len(x))), " i*ny+j identity?", np.array_equal(i * len(ys) + j, np.arange(len(x))))
print("params N_A_ion n_i mu_n mu_p taun taup V_t eps", d["params"])
