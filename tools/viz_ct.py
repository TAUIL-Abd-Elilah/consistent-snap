"""CT cross-section with winding polylines: one z slice of the public CT (level 2 = fit frame), a crop, and any number of
tifxyz sets drawn in given colours. Sets are '<dir>:<b,g,r>[:<scale>]' (scale divides the set's coordinates, e.g. 4 for
L0 human wraps in an L2 frame).
usage: viz_ct.py <z> <y0> <y1> <x0> <x1> <out.png> <set> [<set> ...] [--upscale 2]"""
import os, sys, argparse
import cv2, numpy as np, tifffile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from consistent_snap import slice_poly
os.environ.setdefault('ROUTE24_CACHE', 'ctcache'); os.environ.setdefault('ROUTE24_CAP_GB', '8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'vendor'))
import zfetch

VOL = 'PHerc0139/volumes/20260102150214-2.399um-0.2m-78keV-masked.zarr'
ap = argparse.ArgumentParser()
ap.add_argument('z', type=int); ap.add_argument('y0', type=int); ap.add_argument('y1', type=int)
ap.add_argument('x0', type=int); ap.add_argument('x1', type=int); ap.add_argument('out')
ap.add_argument('sets', nargs='+'); ap.add_argument('--upscale', type=int, default=2); ap.add_argument('--level', type=int, default=2)
ap.add_argument('--label', default='')
a = ap.parse_args()
ct = zfetch.read_box(VOL, a.level, (a.z, a.y0, a.x0), (a.z + 1, a.y1, a.x1))[0].astype(np.float32)
lo, hi = np.percentile(ct[ct > 0], [1, 99.5]) if (ct > 0).any() else (0, 255)
ct = np.clip((ct - lo) / max(hi - lo, 1) * 255, 0, 255).astype(np.uint8)
u = a.upscale
img = cv2.cvtColor(cv2.resize(ct, (ct.shape[1] * u, ct.shape[0] * u), interpolation=cv2.INTER_LINEAR), cv2.COLOR_GRAY2BGR)
panels = [img.copy()]
for spec in a.sets:
    parts = spec.split(':')
    if len(parts) >= 3 and ',' in parts[-2]:
        d, col, sc = ':'.join(parts[:-2]), parts[-2], float(parts[-1])
    else:
        d, col, sc = ':'.join(parts[:-1]), parts[-1], 1.0
    col = tuple(int(c) for c in col.split(','))
    pan = img.copy()
    for name in sorted(os.listdir(d)):
        p = os.path.join(d, name)
        if 'spliced' in name or not os.path.exists(os.path.join(p, 'x.tif')):
            continue
        P = np.stack([tifffile.imread(os.path.join(p, f'{c}.tif')).astype(np.float32) for c in 'xyz'], -1)
        P[~(P > -0.5).all(-1)] = np.nan
        P = P / sc
        L = slice_poly(P, float(a.z)) - [a.x0, a.y0]
        for q0, q1 in zip(L[:-1], L[1:]):
            if np.isfinite(q0).all() and np.isfinite(q1).all() and np.hypot(*(q1 - q0)) < 40:
                cv2.line(pan, tuple(np.rint(q0 * u).astype(int)), tuple(np.rint(q1 * u).astype(int)), col, 1, cv2.LINE_AA)
    panels.append(pan)
out = np.concatenate(panels, axis=1)
if a.label:
    cv2.putText(out, a.label, (8, 22), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 255, 255), 1, cv2.LINE_AA)
cv2.imwrite(a.out, out)
print(a.out, out.shape)
