"""Split the gen-19 full-scroll mesh (one tifxyz, all windings) into one tifxyz per winding using windings.tif, keeping
only rows that reach the band. Ambiguous (-2) and invalid (-1) nodes are dropped. Also checks that the labels belong to
this grid orientation: median radius from the umbilicus must rise with the winding id.
usage: split_fit.py <dl dir> <umbilicus.json> <out dir> z0 z1"""
import json, os, sys
import numpy as np, tifffile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from consistent_snap import Umbilicus

d, umb_p, out, z0, z1 = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5])
X, Y, Z = (tifffile.imread(os.path.join(d, f'{c}.tif')) for c in 'xyz')
W = tifffile.imread(os.path.join(d, 'windings.tif'))
print('grid', X.shape, W.shape, flush=True)
assert X.shape == W.shape
valid = (X > -0.5) & (Y > -0.5) & (Z > -0.5)
umb = Umbilicus(umb_p)
inband = valid & (Z >= z0) & (Z < z1)
lab = np.where(inband, W, -1).astype(np.int32)
# orientation check on the band
rows = []
for w in range(10, 136):
    m = lab == w
    if m.sum() < 500:
        continue
    zz = Z[m]
    cx = np.interp(zz, umb.z, umb.x); cy = np.interp(zz, umb.z, umb.y)
    r = np.hypot(X[m] - cx, Y[m] - cy)
    rows.append((w, int(m.sum()), float(np.median(r))))
ws = np.array([r[0] for r in rows]); rs = np.array([r[2] for r in rows])
rho = np.corrcoef(np.argsort(np.argsort(ws)), np.argsort(np.argsort(rs)))[0, 1]
print('windings in band', len(rows), 'spearman(w, median radius)', round(float(rho), 4))
print('first/last', rows[:3], rows[-3:])
mono = float(np.mean(np.diff(rs) > 0))
print('share of consecutive windings with larger median radius', round(mono, 3))
if rho < 0.95:
    sys.exit('labels do not follow radius: orientation mismatch, stop')
os.makedirs(out, exist_ok=True)
meta = {'type': 'seg', 'format': 'tifxyz', 'scale': [0.05, 0.05], 'source': 'PHercParis4 spiral-output-generation-19 split by windings.tif'}
n_nodes = 0
for w, n, _ in rows:
    m = lab == w
    rr = np.nonzero(m.any(1))[0]; cc = np.nonzero(m.any(0))[0]
    r0, r1, c0, c1 = rr.min(), rr.max() + 1, cc.min(), cc.max() + 1
    # keep the label's rows/cols plus nodes of this winding just outside the band (one row margin) so slices at z0 work
    mm = (W[r0:r1, c0:c1] == w) & valid[r0:r1, c0:c1]
    sub = [np.where(mm, A[r0:r1, c0:c1], -1.0).astype(np.float32) for A in (X, Y, Z)]
    p = os.path.join(out, f'w{w:03d}_g19')
    os.makedirs(p, exist_ok=True)
    for c, A in zip('xyz', sub):
        tifffile.imwrite(os.path.join(p, f'{c}.tif'), A)
    json.dump(dict(meta, grid_rows=[int(r0), int(r1)], grid_cols=[int(c0), int(c1)]), open(os.path.join(p, 'meta.json'), 'w'))
    n_nodes += int(mm.sum())
print('wrote', len(rows), 'windings,', n_nodes, 'nodes')
