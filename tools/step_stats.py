"""Truth-free step count on a snapped run: share of in-slice neighbour pairs (grid axis 1) and z neighbour pairs
(axis 0), both moved/valid, whose snap displacements differ by more than half the local winding pitch.
usage: step_stats.py <snapped dir> [<snapped dir> ...]"""
import os, sys, json, re
import numpy as np, tifffile
zr = None
args = sys.argv[1:]
if args and args[0] == '--z':
    zr = (float(args[1]), float(args[2])); args = args[3:]
for d in args:
    W = {}
    for name in sorted(os.listdir(d)):
        p = os.path.join(d, name)
        m = re.match(r'w(\d+)', name)
        if not m or not os.path.exists(os.path.join(p, 'snap_d.tif')):
            continue
        W[int(m.group(1))] = (tifffile.imread(os.path.join(p, 'snap_d.tif')).astype(np.float32),
                              tifffile.imread(os.path.join(p, 'snap_flag.tif')),
                              tifffile.imread(os.path.join(p, 'x.tif')) > -0.5, tifffile.imread(os.path.join(p, 'z.tif')))
    tot = {'ax1_pairs': 0, 'ax1_steps': 0, 'ax0_pairs': 0, 'ax0_steps': 0}
    for w, (D, F, V, Z) in W.items():
        moved = V & (F != 3) & np.isfinite(D)
        if zr is not None:
            moved &= (Z >= zr[0]) & (Z < zr[1])
        for ax, key in ((1, 'ax1'), (0, 'ax0')):
            a = moved[:, :-1] & moved[:, 1:] if ax == 1 else moved[:-1, :] & moved[1:, :]
            dd = np.abs(np.diff(D, axis=ax))[a]
            tot[f'{key}_pairs'] += int(a.sum()); tot[f'{key}_steps'] += 0
            tot.setdefault(f'{key}_dd', []).append(dd)
    out = {}
    for key in ('ax1', 'ax0'):
        dd = np.concatenate(tot[f'{key}_dd'])
        out[key] = {'pairs': int(len(dd)), 'gt3': float((dd > 3).mean()), 'gt5': float((dd > 5).mean()),
                    'gt8': float((dd > 8).mean()), 'p50': float(np.median(dd)), 'p99': float(np.percentile(dd, 99))}
    print(d, json.dumps(out))
