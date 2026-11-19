"""Duplicate-sheet share exactly as the November scorer computes it (score_ab.duplicate_sheet_share, unchanged), for
any meshes directories. usage: dup_native.py NAME=<dir> ..."""
import os, sys, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'vendor'))
import score_ab as sa
for spec in sys.argv[1:]:
    name, d = spec.split('=', 1)
    t0 = time.time()
    r = sa.duplicate_sheet_share(sa.load_windings(d), 9.596)
    print(name, json.dumps({k: v for k, v in r.items() if not isinstance(v, (list, dict))}), f'{time.time() - t0:.0f} s', flush=True)
