"""Draw one z slice: surface prediction (grey) + fitted windings (red) + snapped windings (green) in a crop.
usage: viz_slice.py <surface zarr array> <z> <out.png> <y0 y1 x0 x1> <meshes dir>:<bgr> [<meshes dir>:<bgr> ...]"""
import os, sys
import cv2, numpy as np, tifffile, zarr
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from consistent_snap import slice_poly

S = zarr.open(sys.argv[1], mode='r'); z = int(sys.argv[2]); out = sys.argv[3]
y0, y1, x0, x1 = map(int, sys.argv[4:8])
img = np.asarray(S[z, y0:y1, x0:x1]).astype(np.uint8)
img = cv2.cvtColor((img // 3).astype(np.uint8), cv2.COLOR_GRAY2BGR)
for spec in sys.argv[8:]:
    d, col = spec.rsplit(':', 1); col = tuple(int(c) for c in col.split(','))
    for name in sorted(os.listdir(d)):
        p = os.path.join(d, name)
        if 'spliced' in name or not os.path.exists(os.path.join(p, 'x.tif')):
            continue
        P = np.stack([tifffile.imread(os.path.join(p, f'{c}.tif')).astype(np.float32) for c in 'xyz'], -1)
        P[~(P > -0.5).all(-1)] = np.nan
        L = slice_poly(P, float(z))
        pts = L - [x0, y0]
        for a, b in zip(pts[:-1], pts[1:]):
            if np.isfinite(a).all() and np.isfinite(b).all() and np.hypot(*(b - a)) < 60:
                cv2.line(img, tuple(np.rint(a).astype(int)), tuple(np.rint(b).astype(int)), col, 1, cv2.LINE_AA)
cv2.imwrite(out, img)
print(out, img.shape)
