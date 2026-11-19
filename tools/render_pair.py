"""Paired surface renders: raw fit (A) vs matched-only snap (Cm) on the same grid pixels, CT level 1 sampled along the
surface normal; human wraps in the same box as a reference. Rules: NOTES.md v4 (fixed before any render).
usage: render_pair.py --box z0 z1 y0 y1 x0 x1 --raw <fit meshes> --cm <Cm meshes> --segs <wrap dir> --out <dir>"""
import argparse, json, os, sys
import cv2, numpy as np, tifffile
from scipy.ndimage import map_coordinates, laplace
os.environ.setdefault('ROUTE24_CACHE', 'ctcache'); os.environ.setdefault('ROUTE24_CAP_GB', '12')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'vendor'))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import zfetch
from consistent_snap import densify

VOL = 'PHerc0139/volumes/20260102150214-2.399um-0.2m-78keV-masked.zarr'
ap = argparse.ArgumentParser()
ap.add_argument('--box', type=int, nargs=6, required=True)
ap.add_argument('--raw', required=True)
ap.add_argument('--cm', required=True)
ap.add_argument('--segs', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--min-nodes', type=int, default=2000)
ap.add_argument('--K', type=int, default=3)
a = ap.parse_args()
z0, z1, y0, y1, x0, x1 = a.box
os.makedirs(a.out, exist_ok=True)
M = 2 * a.K + 4  # margin in L1 voxels
lo = np.array([2 * z0 - M, 2 * y0 - M, 2 * x0 - M]); hi = np.array([2 * z1 + M, 2 * y1 + M, 2 * x1 + M])
ct = zfetch.read_box(VOL, 1, tuple(lo), tuple(hi))
print('ct L1 box', ct.shape, flush=True)
vlo, vhi = np.percentile(ct[::4, ::4, ::4], [1, 99])
med = float(np.median(ct[::4, ::4, ::4]))
offs = np.arange(-a.K, a.K + 1, dtype=np.float32)


def load(p, div=1.0):
    G = np.stack([tifffile.imread(os.path.join(p, f'{c}.tif')).astype(np.float32) for c in 'zyx'], -1)
    G[~(G > -0.5).all(-1)] = np.nan
    return G / div


def inbox(G):
    return np.isfinite(G).all(-1) & (G[..., 0] >= z0) & (G[..., 0] < z1) & (G[..., 1] >= y0) & (G[..., 1] < y1) \
        & (G[..., 2] >= x0) & (G[..., 2] < x1)


def upsample(G, f):
    return densify(G, f)


def render(G):
    """G: grid in L2 coords (NaN = invalid). Returns stack I[off, H, W] (NaN outside)."""
    du = np.gradient(G, axis=0); dv = np.gradient(G, axis=1)
    n = np.cross(du, dv); ln = np.linalg.norm(n, axis=-1, keepdims=True)
    n = n / np.where(ln > 1e-9, ln, np.nan)
    ok = np.isfinite(G).all(-1) & np.isfinite(n).all(-1)
    P = G[ok] * 2 - lo  # L1 box coords
    N = n[ok]
    S = P[None] + offs[:, None, None] * N[None]
    I = map_coordinates(ct, S.reshape(-1, 3).T, order=1, mode='nearest').reshape(len(offs), -1)
    out = np.full((len(offs),) + G.shape[:2], np.nan, np.float32)
    out[:, ok] = I
    return out


def stats(img, m):
    v = img[m]
    lap = laplace(np.nan_to_num(img, nan=med).astype(np.float32))
    er = cv2.erode(m.astype(np.uint8), np.ones((3, 3), np.uint8)).astype(bool)
    return dict(mean0=float(v.mean()), above_median=float((v > med).mean()), texture=float(lap[er].std()), n=int(m.sum()))


def to8(img, m):
    o = np.clip((np.nan_to_num(img, nan=vlo) - vlo) / (vhi - vlo) * 255, 0, 255).astype(np.uint8)
    o[~m] = 0
    return o


res = {'box': a.box, 'ct_median': med, 'window': [float(vlo), float(vhi)], 'windings': {}}
pool = {'A': [], 'Cm': []}
for wn in sorted(os.listdir(a.cm)):
    pc, pr = os.path.join(a.cm, wn), os.path.join(a.raw, wn)
    if not os.path.exists(os.path.join(pc, 'x.tif')) or not os.path.exists(os.path.join(pr, 'x.tif')):
        continue
    C = load(pc)
    ib = inbox(C)
    if ib.sum() < a.min_nodes:
        continue
    A = densify(load(pr), 2)
    if A.shape != C.shape:
        print('shape mismatch', wn, A.shape, C.shape); continue
    r, c = np.where(ib)
    r0, r1, c0, c1 = r.min(), r.max() + 1, c.min(), c.max() + 1
    Cs, As = C[r0:r1, c0:c1].copy(), A[r0:r1, c0:c1].copy()
    Cs[~ib[r0:r1, c0:c1]] = np.nan
    step = np.nanmedian(np.linalg.norm(np.diff(As, axis=1), axis=-1))
    f = max(1, int(round(step * 2)))  # ~1 L1 voxel per pixel
    Cu, Au = upsample(Cs, f), upsample(As, f)
    m = inbox(Cu) & np.isfinite(Au).all(-1)
    if m.sum() < 1000:
        continue
    IA, IC = render(Au), render(Cu)
    m &= np.isfinite(IA[a.K]) & np.isfinite(IC[a.K])
    sA, sC = stats(IA[a.K], m), stats(IC[a.K], m)
    res['windings'][wn] = dict(nodes=int(ib.sum()), f=f, A=sA, Cm=sC)
    pool['A'].append(IA[a.K][m]); pool['Cm'].append(IC[a.K][m])
    print(wn, 'nodes', int(ib.sum()), 'px', int(m.sum()), 'A', {k: round(v, 3) for k, v in sA.items()},
          'Cm', {k: round(v, 3) for k, v in sC.items()}, flush=True)
    for name, I in (('A', IA), ('Cm', IC)):
        cv2.imwrite(os.path.join(a.out, f'{wn}_{name}_0.png'), to8(I[a.K], m))
        cv2.imwrite(os.path.join(a.out, f'{wn}_{name}_max.png'), to8(np.nanmax(I, 0), m))
for k in pool:
    v = np.concatenate(pool[k]) if pool[k] else np.zeros(1)
    res['pooled_' + k] = dict(mean0=float(v.mean()), above_median=float((v > med).mean()), n=int(v.size))
ws = res['windings']
res['A_ge_Cm'] = {s: [w for w in ws if ws[w]['A'][s] >= ws[w]['Cm'][s]] for s in ('mean0', 'above_median', 'texture')}
res['top3'] = sorted(ws, key=lambda w: -ws[w]['Cm']['n'])[:3]

# human wraps in the same box (reference, own parametrization)
res['wraps'] = {}
for sd in sorted(os.listdir(a.segs)):
    p = os.path.join(a.segs, sd)
    if not os.path.exists(os.path.join(p, 'x.tif')):
        continue
    G = load(p, 4.0)
    ib = inbox(G)
    if ib.sum() < 500:
        continue
    r, c = np.where(ib)
    Gs = G[r.min():r.max() + 1, c.min():c.max() + 1].copy()
    Gs[~ib[r.min():r.max() + 1, c.min():c.max() + 1]] = np.nan
    step = np.nanmedian(np.linalg.norm(np.diff(Gs, axis=1), axis=-1))
    Gu = upsample(Gs, max(1, int(round(step * 2))))
    m = inbox(Gu)
    I = render(Gu)
    m &= np.isfinite(I[a.K])
    tag = sd.split('-')[-1].split('_')[0]
    res['wraps'][tag] = stats(I[a.K], m)
    cv2.imwrite(os.path.join(a.out, f'wrap_{tag}_0.png'), to8(I[a.K], m))
    cv2.imwrite(os.path.join(a.out, f'wrap_{tag}_max.png'), to8(np.nanmax(I, 0), m))
    print('wrap', tag, res['wraps'][tag], flush=True)
json.dump(res, open(os.path.join(a.out, 'render_stats.json'), 'w'), indent=1)
print('pooled A', res['pooled_A'], 'Cm', res['pooled_Cm'])
print('A >= Cm:', res['A_ge_Cm'], 'top3', res['top3'])
