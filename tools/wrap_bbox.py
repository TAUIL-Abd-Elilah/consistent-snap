"""Bounding box (L2 frame) of the human wrap nodes near a given z. usage: wrap_bbox.py <segments dir> <z> [<z> ...]"""
import os, sys
import numpy as np, tifffile

d = sys.argv[1]
zs = [float(z) for z in sys.argv[2:]]
pts = []
for s in sorted(os.listdir(d)):
    p = os.path.join(d, s)
    if not os.path.exists(os.path.join(p, 'x.tif')):
        continue
    P = np.stack([tifffile.imread(os.path.join(p, f'{c}.tif')).astype(np.float32) for c in 'xyz'], -1) / 4.0
    v = (P > -0.1).all(-1)
    pts.append(P[v])
P = np.concatenate(pts)
for z in zs:
    q = P[np.abs(P[:, 2] - z) < 1.0]
    print(z, len(q), 'x', np.percentile(q[:, 0], [1, 50, 99]).round(), 'y', np.percentile(q[:, 1], [1, 50, 99]).round())
