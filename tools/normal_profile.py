"""CT intensity profile across a surface: for valid nodes of every winding inside a box, sample the public CT (level 2,
the fit frame) along the node normal at offsets -R..R voxels. A surface that sits on a sheet shows a peak at 0 that falls
off on both sides ("one voxel off a perfectly snapped surface should lie in air", sean, 3 Nov). Reports the mean profile,
the peak-to-flank contrast I(0) - mean(I(+-4..+-6)), and the share of nodes whose own intensity is the local maximum
within +-2 voxels.
usage: normal_profile.py --box z0 z1 y0 y1 x0 x1 [--max 200000] [--R 8] NAME=<meshes dir>[:scale] ..."""
import argparse, json, os, re, sys
import numpy as np, tifffile
from scipy.ndimage import map_coordinates
os.environ.setdefault('ROUTE24_CACHE', 'ctcache'); os.environ.setdefault('ROUTE24_CAP_GB', '8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'vendor'))
import zfetch

VOL = 'PHerc0139/volumes/20260102150214-2.399um-0.2m-78keV-masked.zarr'
ap = argparse.ArgumentParser()
ap.add_argument('--box', type=int, nargs=6, required=True)
ap.add_argument('--max', type=int, default=200000)
ap.add_argument('--R', type=int, default=8)
ap.add_argument('--vol', default=VOL)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('arms', nargs='+')
a = ap.parse_args()
z0, z1, y0, y1, x0, x1 = a.box
R = a.R
ct = zfetch.read_box(a.vol, 2, (z0 - R - 2, y0 - R - 2, x0 - R - 2), (z1 + R + 2, y1 + R + 2, x1 + R + 2)).astype(np.float32)
org = np.array([z0 - R - 2, y0 - R - 2, x0 - R - 2], np.float32)
print('ct box', ct.shape, 'mean', float(ct.mean()), flush=True)
offs = np.arange(-R, R + 1, dtype=np.float32)
rng = np.random.default_rng(a.seed)
out = {}
for spec in a.arms:
    name, rest = spec.split('=', 1)
    m = re.match(r'^(.*?)(?::([0-9.]+))?$', rest)
    d, sc = m.group(1), float(m.group(2) or 1.0)
    P, N = [], []
    for wn in sorted(os.listdir(d)):
        p = os.path.join(d, wn)
        if 'spliced' in wn or not os.path.exists(os.path.join(p, 'x.tif')):
            continue
        G = np.stack([tifffile.imread(os.path.join(p, f'{c}.tif')).astype(np.float32) for c in 'zyx'], -1)
        v = (G > -0.5).all(-1) & np.isfinite(G).all(-1)
        G = G / sc
        du, dv = np.gradient(G, axis=0), np.gradient(G, axis=1)
        n = np.cross(du, dv)
        ln = np.linalg.norm(n, axis=-1)
        vv = v.copy()
        vv[1:-1, 1:-1] &= v[:-2, 1:-1] & v[2:, 1:-1] & v[1:-1, :-2] & v[1:-1, 2:]
        vv[[0, -1], :] = False
        vv[:, [0, -1]] = False
        ok = vv & (ln > 1e-6)
        q = G[ok]
        inb = (q[:, 0] >= z0) & (q[:, 0] < z1) & (q[:, 1] >= y0) & (q[:, 1] < y1) & (q[:, 2] >= x0) & (q[:, 2] < x1)
        if inb.any():
            P.append(q[inb]); N.append((n[ok] / ln[ok][:, None])[inb])
    P, N = np.concatenate(P), np.concatenate(N)
    if len(P) > a.max:
        sel = rng.choice(len(P), a.max, replace=False)
        P, N = P[sel], N[sel]
    S = P[:, None, :] + offs[None, :, None] * N[:, None, :] - org
    I = map_coordinates(ct, S.reshape(-1, 3).T, order=1, mode='nearest').reshape(len(P), len(offs))
    I = I.astype(np.float64)
    c0 = R
    prof = I.mean(0)
    flank = np.r_[c0 - 6:c0 - 3, c0 + 4:c0 + 7]
    contrast = I[:, c0] - I[:, flank].mean(1)
    near = I[:, c0 - 2:c0 + 3]
    is_max = near.argmax(1) == 2
    out[name] = dict(nodes=int(len(P)), profile=[round(float(x), 1) for x in prof],
                     contrast_mean=round(float(contrast.mean()), 2), contrast_median=round(float(np.median(contrast)), 2),
                     centre_is_local_max=round(float(is_max.mean()), 3), I0_mean=round(float(I[:, c0].mean()), 1),
                     I_pm2_mean=round(float(np.r_[I[:, c0 - 2], I[:, c0 + 2]].mean()), 1))
    print(name, json.dumps({k: v for k, v in out[name].items() if k != 'profile'}), flush=True)
    print('   profile', out[name]['profile'], flush=True)
