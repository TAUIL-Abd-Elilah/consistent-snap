"""Mirror the compressed chunks of one level of a public open-data zarr v2 array that cover a zyx box into a local
zarr directory with the same .zarray, so zarr can open it locally. Resumable (existing chunk files are skipped).

usage: mirror_chunks.py <s3 key of the zarr (no bucket)> <level> <out zarr dir> z0 z1 y0 y1 x0 x1 [--workers 16]
"""
import argparse, json, os, sys, time, urllib.error, urllib.request
from concurrent.futures import ThreadPoolExecutor

BUCKET = 'https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com'


def get(url, tries=6):
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):
                return None
            if a == tries - 1:
                raise
        except (urllib.error.URLError, ConnectionError, TimeoutError, OSError):
            if a == tries - 1:
                raise
        time.sleep(2.0 * (a + 1))


ap = argparse.ArgumentParser()
ap.add_argument('key'); ap.add_argument('level'); ap.add_argument('out')
ap.add_argument('box', nargs=6, type=int); ap.add_argument('--workers', type=int, default=16)
a = ap.parse_args()
z0, z1, y0, y1, x0, x1 = a.box
os.makedirs(f'{a.out}/{a.level}', exist_ok=True)
for rel in ('.zgroup', '.zattrs'):
    raw = get(f'{BUCKET}/{a.key}/{rel}')
    if raw is not None:
        open(f'{a.out}/{rel}', 'wb').write(raw)
zarray = get(f'{BUCKET}/{a.key}/{a.level}/.zarray')
open(f'{a.out}/{a.level}/.zarray', 'wb').write(zarray)
m = json.loads(zarray)
cz, cy, cx = m['chunks']
keys = [(i, j, k) for i in range(z0 // cz, (z1 - 1) // cz + 1) for j in range(y0 // cy, (y1 - 1) // cy + 1)
        for k in range(x0 // cx, (x1 - 1) // cx + 1)]
sep = m.get('dimension_separator', '.')
done = {'n': 0, 'bytes': 0, 'missing': 0}


def fetch(key):
    rel = sep.join(map(str, key))
    path = os.path.join(a.out, a.level, *rel.split('/'))
    if os.path.exists(path):
        done['n'] += 1
        return
    raw = get(f'{BUCKET}/{a.key}/{a.level}/{rel}')
    if raw is None:                      # absent chunk = fill value
        done['missing'] += 1; done['n'] += 1
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + '.part'
    open(tmp, 'wb').write(raw)
    for t in range(20):                  # Windows: a scanner may hold the new file for a moment
        try:
            os.replace(tmp, path); break
        except PermissionError:
            time.sleep(0.5 * (t + 1))
    done['n'] += 1; done['bytes'] += len(raw)
    if done['n'] % 100 == 0:
        print(f"{done['n']}/{len(keys)} chunks, {done['bytes'] / 1e6:.0f} MB, {done['missing']} absent", flush=True)


print(f'{len(keys)} chunks to mirror', flush=True)
with ThreadPoolExecutor(a.workers) as pool:
    list(pool.map(fetch, keys))
print(f"DONE {done['n']}/{len(keys)} chunks, {done['bytes'] / 1e6:.0f} MB, {done['missing']} absent", flush=True)
