"""Cut the stair-steps out of a snapped run: a node whose radial move differs from a 4-neighbour's by more than --thr
voxels sits on a sheet jump (the connector between two sheets runs through air), so both nodes of every such pair become
invalid (-1). Pieces smaller than --min-piece nodes left between cuts are dropped too. Writes a new run directory and a
report of the area kept. usage: cut_steps.py <snapped dir> <out dir> [--thr 6] [--min-piece 50]"""
import argparse, json, os, re, shutil
import numpy as np, tifffile
from scipy.ndimage import label

ap = argparse.ArgumentParser()
ap.add_argument('src'); ap.add_argument('dst')
ap.add_argument('--thr', type=float, default=6.0)
ap.add_argument('--min-piece', type=int, default=50)
a = ap.parse_args()
os.makedirs(a.dst, exist_ok=True)
rep = {'thr_vx': a.thr, 'min_piece': a.min_piece, 'windings': {}}
tot = np.zeros(3, np.int64)
for name in sorted(os.listdir(a.src)):
    p = os.path.join(a.src, name)
    if not re.match(r'w\d+', name) or not os.path.exists(os.path.join(p, 'x.tif')):
        continue
    X, Y, Z = (tifffile.imread(os.path.join(p, f'{c}.tif')) for c in 'xyz')
    D = tifffile.imread(os.path.join(p, 'snap_d.tif')).astype(np.float32)
    F = tifffile.imread(os.path.join(p, 'snap_flag.tif'))
    V = (X > -0.5) & (Y > -0.5) & (Z > -0.5)
    mv = V & (F != 3)
    bad = np.zeros_like(V)
    for ax in (0, 1):
        sl0 = (slice(None, -1), slice(None)) if ax == 0 else (slice(None), slice(None, -1))
        sl1 = (slice(1, None), slice(None)) if ax == 0 else (slice(None), slice(1, None))
        both = mv[sl0] & mv[sl1]
        jump = both & (np.abs(D[sl0] - D[sl1]) > a.thr)
        bad[sl0] |= jump
        bad[sl1] |= jump
    keep = V & ~bad
    lab, n = label(keep)
    if n:
        sizes = np.bincount(lab.ravel())
        small = sizes < a.min_piece
        small[0] = False
        keep &= ~small[lab]
    q = os.path.join(a.dst, name)
    os.makedirs(q, exist_ok=True)
    for c, arr in zip('xyz', (X, Y, Z)):
        tifffile.imwrite(os.path.join(q, f'{c}.tif'), np.where(keep, arr, -1).astype(np.float32))
    tifffile.imwrite(os.path.join(q, 'snap_d.tif'), np.where(keep, D, 0).astype(np.float32))
    tifffile.imwrite(os.path.join(q, 'snap_flag.tif'), np.where(keep, F, 3).astype(np.uint8))
    mp = os.path.join(p, 'meta.json')
    if os.path.exists(mp):
        m = json.load(open(mp))
        m['cut_steps'] = {'thr_vx': a.thr, 'min_piece': a.min_piece, 'of': p}
        json.dump(m, open(os.path.join(q, 'meta.json'), 'w'), indent=1)
    rep['windings'][name] = dict(valid=int(V.sum()), cut=int((V & bad).sum()), kept=int(keep.sum()))
    tot += [int(V.sum()), int((V & bad).sum()), int(keep.sum())]
rep['total'] = dict(valid=int(tot[0]), cut=int(tot[1]), kept=int(tot[2]), kept_share=round(float(tot[2] / max(tot[0], 1)), 4))
json.dump(rep, open(os.path.join(a.dst, 'cut_report.json'), 'w'), indent=1)
print(a.dst, json.dumps(rep['total']))
