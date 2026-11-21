"""Paired ink renders of the raw gen-19 fit (A) and its matched-only snap (Cm) through the public ink-3d prediction,
emulating spiral-fitting/render_ink.py (vc_render_tifxyz --scale 0.25 --num-slices 5, max composite, /p95 -> uint8).
Rules: PREREG_P4.md. Both arms are rendered on the gen-19 global (flattened) grid, upsampled x5 (1 px = 4 level-2
voxels), on the pixels valid in Cm only. Ink is read chunk by chunk (never the whole band in RAM).
usage: render_ink_pair.py <fit_A dir> <Cm dir> <out dir>"""
import os, sys, json, time
import numpy as np, tifffile, cv2
from scipy.ndimage import map_coordinates
os.environ.setdefault('ROUTE24_CACHE', 'data/cache'); os.environ.setdefault('ROUTE24_CAP_GB', '16')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'vendor'))
import zfetch

INK = 'PHercParis4/representations/predictions/ink-3d/20260411134726-ink3d-20260926084000-student-v12A-25000-c7766be9-L1.zarr'
LEVEL = 1          # = level 2 of the scan (9.6 um), the fit frame
UP = 5             # native grid step 20 vx -> 4 vx per pixel (render_ink --scale 0.25)
OFFS = np.arange(-2, 3, dtype=np.float32)   # --num-slices 5, 1 voxel apart
fa, fc, out = sys.argv[1:4]
INK_LOCAL = os.environ.get('INK_LOCAL', 'data/p4_ink_band.zarr/1')  # tools/mirror_chunks.py of INK, level 1
os.makedirs(out, exist_ok=True)
t0 = time.time()


def load(p):
    G = np.stack([tifffile.imread(os.path.join(p, f'{c}.tif')).astype(np.float32) for c in 'zyx'], -1)
    G[~(G > -0.5).all(-1)] = np.nan
    return G


def resample(G, step, H, W):
    """Sample grid G (node spacing = step native cells) at native coords i/UP, j/UP for i<H, j<W (bilinear, valid
    only where all four parents are valid)."""
    ys = np.arange(H, dtype=np.float32) / UP / step
    xs = np.arange(W, dtype=np.float32) / UP / step
    mx, my = np.meshgrid(xs, ys)
    v = np.isfinite(G).all(-1).astype(np.float32)
    o = np.stack([cv2.remap(np.nan_to_num(G[..., c]), mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT,
                            borderValue=0) for c in range(3)], -1)
    vv = cv2.remap(v, mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT, borderValue=0)
    o[vv < 0.999] = np.nan
    return o


def normals(P):
    du = np.gradient(P, axis=0); dv = np.gradient(P, axis=1)
    n = np.cross(du, dv); ln = np.linalg.norm(n, axis=-1, keepdims=True)
    return n / np.where(ln > 1e-9, ln, np.nan)


# ---- 1. pixel positions and normals for both arms, per winding, on the global grid
r_min, r_max, c_min, c_max = 10 ** 9, 0, 10 ** 9, 0
wins = []
for name in sorted(os.listdir(fc)):
    pc = os.path.join(fc, name)
    if not os.path.exists(os.path.join(pc, 'x.tif')):
        continue
    meta = json.load(open(os.path.join(fa, name, 'meta.json')))
    r0, r1 = meta['grid_rows']; c0, c1 = meta['grid_cols']
    wins.append((name, r0, r1, c0, c1))
    r_min, r_max, c_min, c_max = min(r_min, r0), max(r_max, r1), min(c_min, c0), max(c_max, c1)
Hg, Wg = (r_max - r_min - 1) * UP + 1, (c_max - c_min - 1) * UP + 1
print(f'global strip {Hg} x {Wg} px, {len(wins)} windings', flush=True)


def build(arm):
    G, Pl, Nl = [], [], []
    for name, r0, r1, c0, c1 in wins:
        A = load(os.path.join(fa, name)); C = load(os.path.join(fc, name))
        H, W = (A.shape[0] - 1) * UP + 1, (A.shape[1] - 1) * UP + 1
        Au, Cu = resample(A, 1, H, W), resample(C, 0.5, H, W)
        m = np.isfinite(Cu).all(-1) & np.isfinite(Au).all(-1)
        nA, nC = normals(Au), normals(Cu)
        m &= np.isfinite(nA).all(-1) & np.isfinite(nC).all(-1)
        ii, jj = np.nonzero(m)
        G.append((ii + (r0 - r_min) * UP).astype(np.int64) * Wg + (jj + (c0 - c_min) * UP))
        Pl.append((Au if arm == 'A' else Cu)[m]); Nl.append((nA if arm == 'A' else nC)[m])
        del A, C, Au, Cu, nA, nC
    return np.concatenate(G), np.concatenate(Pl), np.concatenate(Nl)


# ---- 2. sample the ink chunk by chunk
import zarr
ZI = zarr.open(INK_LOCAL, mode='r')
ch = np.asarray(ZI.chunks)
shp = np.asarray(ZI.shape)
img = {}
for arm in ('A', 'Cm'):
    gi, Pall, Nall = build(arm)
    print(arm, 'pixels', len(gi), f'{time.time() - t0:.0f} s', flush=True)
    best = np.zeros(len(gi), np.float32)
    for b0 in range(0, len(gi), 12_000_000):
      P = Pall[b0:b0 + 12_000_000]; N = Nall[b0:b0 + 12_000_000]
      bbest = np.zeros(len(P), np.float32)
      for o in OFFS:
          S = P + o * N
          key = np.floor(S / ch).astype(np.int32)
          kid = (key[:, 0].astype(np.int64) * 4096 + key[:, 1]) * 4096 + key[:, 2]
          del S
          S = None
          order = np.argsort(kid, kind='stable'); ks = kid[order]
          starts = np.r_[0, np.nonzero(np.diff(ks))[0] + 1]; ends = np.r_[starts[1:], len(ks)]
          val = np.zeros(len(P), np.float32)
          for s0, s1 in zip(starts, ends):
              idx = order[s0:s1]
              k = key[idx[0]]
              lo = k * ch - 2; hi = (k + 1) * ch + 2
              lo_c = np.maximum(lo, 0); hi_c = np.minimum(hi, shp)
              box = ZI[lo_c[0]:hi_c[0], lo_c[1]:hi_c[1], lo_c[2]:hi_c[2]].astype(np.float32)
              val[idx] = map_coordinates(box, (P[idx] + o * N[idx] - lo_c).T, order=1, mode='nearest')
          bbest = np.maximum(bbest, val)
          print(arm, 'batch', b0, 'offset', o, f'{len(starts)} chunks, {time.time() - t0:.0f} s', flush=True)
      best[b0:b0 + len(P)] = bbest
    im = np.zeros(Hg * Wg, np.float32); im[gi] = best
    img[arm] = im.reshape(Hg, Wg)
    del Pall, Nall, best, im
    np.save(os.path.join(out, f'{arm}_raw.npy'), img[arm].astype(np.float16))
paired = np.zeros(Hg * Wg, bool); paired[gi] = True; paired = paired.reshape(Hg, Wg)
np.save(os.path.join(out, 'paired_mask.npy'), paired)
stats = {}
for arm in ('A', 'Cm'):
    s = img[arm]
    p95 = float(np.percentile(s, 95))
    s8 = (np.clip(s / p95, 0, 1) * 255).astype(np.uint8) if p95 > 0 else s.astype(np.uint8)
    d = os.path.join(out, arm, 'ink'); os.makedirs(d, exist_ok=True)
    for t in range(0, Wg, 16384):
        cv2.imwrite(os.path.join(d, f'w010-135_flat.{t // 16384:03d}.jpg'), s8[:, t:t + 16384], [cv2.IMWRITE_JPEG_QUALITY, 95])
    stats[arm] = dict(p95=p95, mean_paired=float(s[paired].mean()), paired_px=int(paired.sum()))
stats['strip'] = [Hg, Wg]
json.dump(stats, open(os.path.join(out, 'render_stats.json'), 'w'), indent=1)
print(stats, f'{time.time() - t0:.0f} s')
