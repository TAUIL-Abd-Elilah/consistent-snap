"""Precision of a surface against the human wraps: among nodes of the arm inside the scored z window that lie within
--zone vx of a human-wrap triangle, the share within 2 vx (and 4 vx). Human wraps: L0 tifxyz / 4 (fit frame), every
--stride-th grid node. Uses score_band.Mesh (point-to-triangle distance, unchanged).
usage: precision.py [--segs DIR] [--z0 3200 --z1 4200] [--max 400000] NAME=<meshes dir> ..."""
import argparse, os, sys, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'vendor'))
import score_ab as sa
from score_band import Mesh
ap = argparse.ArgumentParser()
ap.add_argument('--segs', required=True, help='human wrap segments (L0 tifxyz), e.g. segments_masked')
ap.add_argument('--z0', type=float, default=3200); ap.add_argument('--z1', type=float, default=4200)
ap.add_argument('--edge', type=float, default=20); ap.add_argument('--zone', type=float, default=12)
ap.add_argument('--stride', type=int, default=1); ap.add_argument('--max', type=int, default=400000)
ap.add_argument('arms', nargs='+')
a = ap.parse_args()
W = []
for i, n in enumerate(sorted(os.listdir(a.segs))):
    p = os.path.join(a.segs, n)
    if not os.path.exists(os.path.join(p, 'x.tif')):
        continue
    zyx, v = sa.load_tifxyz(p)
    zyx = np.where(v[..., None], zyx / 4.0, -1.0)[::a.stride, ::a.stride]
    W.append((i, n, zyx, v[::a.stride, ::a.stride]))
mesh = Mesh(W, max_dist=a.zone)
print('wrap triangles', len(mesh.tris), 'dropped quads', mesh.dropped, flush=True)
rng = np.random.default_rng(0)
for spec in a.arms:
    name, d = spec.split('=', 1)
    P = []
    for w, n, zyx, v in sa.load_windings(d):
        q = zyx[v]
        P.append(q[(q[:, 0] >= a.z0 + a.edge) & (q[:, 0] < a.z1 - a.edge)])
    P = np.concatenate(P)
    if len(P) > a.max:
        P = P[rng.choice(len(P), a.max, replace=False)]
    dist, _ = mesh.query(P)
    zone = np.isfinite(dist)
    print(name, json.dumps(dict(nodes=int(len(P)), in_zone=int(zone.sum()),
                                prec2=round(float((dist[zone] <= 2).mean()), 4), prec4=round(float((dist[zone] <= 4).mean()), 4))), flush=True)
