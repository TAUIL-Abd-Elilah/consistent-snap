"""Take every f-th node of each snapped winding (corner-aligned densify -> the fit's own nodes after the move), so
density-sensitive metrics compare like with like. usage: native_spacing.py <snapped dir> <out dir> [f=2]"""
import os, sys, json, shutil
import tifffile
src, dst = sys.argv[1], sys.argv[2]
f = int(sys.argv[3]) if len(sys.argv) > 3 else 2
os.makedirs(dst, exist_ok=True)
n = 0
for name in sorted(os.listdir(src)):
    p = os.path.join(src, name)
    if not os.path.isdir(p) or not os.path.exists(os.path.join(p, 'x.tif')):
        continue
    q = os.path.join(dst, name); os.makedirs(q, exist_ok=True)
    for c in ('x', 'y', 'z', 'snap_flag', 'snap_d'):
        fp = os.path.join(p, f'{c}.tif')
        if os.path.exists(fp):
            tifffile.imwrite(os.path.join(q, f'{c}.tif'), tifffile.imread(fp)[::f, ::f])
    mp = os.path.join(p, 'meta.json')
    if os.path.exists(mp):
        m = json.load(open(mp))
        if 'scale' in m:
            m['scale'] = [s / f for s in m['scale']] if isinstance(m['scale'], list) else m['scale'] / f
        m['native_spacing_of'] = p
        json.dump(m, open(os.path.join(q, 'meta.json'), 'w'), indent=1)
    n += 1
print(dst, n, 'windings')
