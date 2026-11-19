"""Keep only nodes with snap_flag in the given set (others -1). usage: mask_flag.py <src> <dst> <flags e.g. 0 or 0,1>"""
import os, sys, re, json
import numpy as np, tifffile
src, dst, keepf = sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(',')]
os.makedirs(dst, exist_ok=True)
kept = tot = 0
for name in sorted(os.listdir(src)):
    p = os.path.join(src, name)
    if not re.match(r'w\d+', name) or not os.path.exists(os.path.join(p, 'x.tif')):
        continue
    F = tifffile.imread(os.path.join(p, 'snap_flag.tif'))
    q = os.path.join(dst, name); os.makedirs(q, exist_ok=True)
    X = tifffile.imread(os.path.join(p, 'x.tif'))
    keep = np.isin(F, keepf) & (X > -0.5)
    kept += int(keep.sum()); tot += int((X > -0.5).sum())
    for c in 'xyz':
        A = tifffile.imread(os.path.join(p, f'{c}.tif'))
        tifffile.imwrite(os.path.join(q, f'{c}.tif'), np.where(keep, A, -1).astype(np.float32))
    for c in ('snap_flag', 'snap_d'):
        tifffile.imwrite(os.path.join(q, f'{c}.tif'), tifffile.imread(os.path.join(p, f'{c}.tif')))
print(dst, 'kept', kept, 'of', tot)
