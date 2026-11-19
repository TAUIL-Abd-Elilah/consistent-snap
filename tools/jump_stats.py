"""Truth-free jump count: in-slice (axis 1) and z (axis 0) neighbour triples of moved nodes; a jump is a second
difference |d[i-1] - 2 d[i] + d[i+1]| > frac * local pitch proxy (here a fixed vx threshold), split by node flags.
usage: jump_stats.py [--z z0 z1] [--thr 4] <snapped dir> ..."""
import os, sys, re, json
import numpy as np, tifffile
args = sys.argv[1:]
zr, thr = None, 4.0
while args and args[0].startswith('--'):
    if args[0] == '--z':
        zr = (float(args[1]), float(args[2])); args = args[3:]
    elif args[0] == '--thr':
        thr = float(args[1]); args = args[2:]
for d in args:
    res = {}
    for ax in (1, 0):
        n_tr = 0; n_j = 0; byflag = {}
        for name in sorted(os.listdir(d)):
            p = os.path.join(d, name)
            if not re.match(r'w\d+', name) or not os.path.exists(os.path.join(p, 'snap_d.tif')):
                continue
            D = tifffile.imread(os.path.join(p, 'snap_d.tif')).astype(np.float32)
            F = tifffile.imread(os.path.join(p, 'snap_flag.tif'))
            V = tifffile.imread(os.path.join(p, 'x.tif')) > -0.5
            Z = tifffile.imread(os.path.join(p, 'z.tif'))
            mv = V & (F != 3)
            if zr is not None:
                mv &= (Z >= zr[0]) & (Z < zr[1])
            if ax == 0:
                D, F, mv = D.T, F.T, mv.T
            a, b, c = mv[:, :-2], mv[:, 1:-1], mv[:, 2:]
            ok = a & b & c
            dd = np.abs(D[:, :-2] - 2 * D[:, 1:-1] + D[:, 2:])[ok]
            fl = F[:, 1:-1][ok]
            n_tr += int(ok.sum()); n_j += int((dd > thr).sum())
            for f in (0, 1):
                byflag.setdefault(f, [0, 0])
                byflag[f][0] += int((fl == f).sum()); byflag[f][1] += int(((fl == f) & (dd > thr)).sum())
        res['inslice' if ax == 1 else 'z'] = dict(triples=n_tr, jump_share=round(n_j / max(n_tr, 1), 4),
                                                  matched_centre=round(byflag[0][1] / max(byflag[0][0], 1), 4),
                                                  filled_centre=round(byflag[1][1] / max(byflag[1][0], 1), 4))
    print(d, json.dumps(res))
