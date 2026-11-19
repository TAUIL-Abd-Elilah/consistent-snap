"""Truth-free diagnostics of a label_snap dump: for the largest components, how the fit's phi drifts against the
sheet-graph label (dev = phi - label - base) over radius and angle. usage: label_diag.py <dump.npz> [top]"""
import sys
import numpy as np

d = np.load(sys.argv[1])
top = int(sys.argv[2]) if len(sys.argv) > 2 else 6
R, PHI, RAY, root, lab, base = d['R'], d['PHI'], d['RAY'], d['root'], d['lab'], d['base']
M = int(d['M'])
dev = PHI - lab - base[root]
sizes = np.bincount(root)
for c in np.argsort(sizes)[::-1][:top]:
    s = root == c
    m = RAY[s] % M
    k = RAY[s] // M
    dv = dev[s]
    f = np.isfinite(dv)
    print(f'comp {c}: {s.sum()} nodes, slices {k.min()}-{k.max()}, rays covered {len(np.unique(m))}, '
          f'radius {np.percentile(R[s], [1, 50, 99]).round()}, label range {lab[s].min()}..{lab[s].max()}')
    print('   dev pct 1/10/50/90/99:', np.percentile(dv[f], [1, 10, 50, 90, 99]).round(2))
    # dev vs label (turn): median dev per relative label
    for L in np.unique(lab[s])[:12]:
        q = s & (lab == L) & np.isfinite(dev)
        if q.sum() > 200:
            print(f'     label {L:+d}: n {q.sum():7d} radius {np.median(R[q]):7.1f} dev median {np.median(dev[q]):+.2f} '
                  f'p10 {np.percentile(dev[q], 10):+.2f} p90 {np.percentile(dev[q], 90):+.2f}')
