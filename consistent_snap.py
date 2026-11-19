"""consistent_snap.py: snap ALL windings of a spiral fit onto a surface prediction at once, so that along every ray from
the scroll axis the windings keep their order and no two windings land on the same predicted sheet.

Why: per-surface snappers move each winding to its nearest predicted sheet independently. Where the fit's pitch does not
match the scroll, two windings snap onto the same sheet (stacking) and a sheet in between is left empty. Here the
assignment is solved jointly per ray with an exact order-preserving alignment (dynamic programming):

    fitted windings along the ray, in radial order, target radii t_1 < t_2 < ... < t_n
    predicted sheets along the ray, in radial order,  rho_1 < rho_2 < ... < rho_m
    cost = sum over matched pairs of ((t_i - rho_j) / (sigma * p_i))^2     (p_i = local winding spacing on this ray)
         + C_w for each winding left unmatched
         + C_s for each sheet left without a winding
    matches are one-to-one and order-preserving, and only allowed when |t_i - rho_j| <= d_max * p_i.

Two passes make the result coherent across neighbouring rays and slices:
  pass 1: t_i = fitted radius r_i  ->  raw radial displacement per winding / slice / ray;
  smooth: per winding, a robust (median, then normalised-Gaussian) displacement field over (slice, ray), gaps filled;
  pass 2: t_i = r_i + smoothed displacement, with a tighter d_max  ->  final matches.
Every node of a winding (grid densified f x) is moved radially by its final displacement: the pass-2 match where its ray
found one (flag 0), otherwise the smoothed field (flag 1); nodes with neither keep their position (flag 3). snap_flag.tif
holds the flag, snap_d.tif the radial move in voxels.

--mode naive replaces the joint alignment by nearest-sheet-per-winding on the same rays and the same two passes
(the per-surface baseline: identical except for the ordering / one-winding-per-sheet constraint).

Inputs: a fit_spiral meshes directory (w###_<tag>/{x,y,z}.tif tifxyz, fit-frame voxels), a surface-prediction zarr array
in the same frame (uint8; value > --thr is sheet), the umbilicus JSON (control points x, y, z in the fit frame).
Output: tifxyz per winding (+ snap_flag.tif, snap_d.tif) and report.json.  CPU only.
"""
from __future__ import annotations

import argparse
import json
import os
import time

import cv2
import numpy as np
import tifffile
import zarr
from numba import njit, prange
from scipy.ndimage import gaussian_filter, median_filter

# --------------------------------------------------------------------------------------------- mesh I/O and geometry


def load_windings(d, include_spliced=False):
    out = {}
    for name in sorted(os.listdir(d)):
        p = os.path.join(d, name)
        if not os.path.exists(os.path.join(p, 'x.tif')) or (('spliced' in name) and not include_spliced):
            continue
        w = int(name[1:4])
        P = np.stack([tifffile.imread(os.path.join(p, f'{c}.tif')).astype(np.float32) for c in 'xyz'], -1)
        P[~(P > -0.5).all(-1)] = np.nan
        mp = os.path.join(p, 'meta.json')
        out[w] = (name, P, json.load(open(mp)) if os.path.exists(mp) else {})
    return out


def densify(P, f):
    """Bilinear upsampling of a tifxyz grid by f (corner-aligned); valid only where all four parent nodes are valid."""
    if f == 1:
        return P
    H, W, _ = P.shape
    v = np.isfinite(P).all(-1).astype(np.float32)
    Pz = np.nan_to_num(P)
    ys = np.arange((H - 1) * f + 1, dtype=np.float32) / f
    xs = np.arange((W - 1) * f + 1, dtype=np.float32) / f
    mx, my = np.meshgrid(xs, ys)
    out = np.stack([cv2.remap(Pz[..., c], mx, my, cv2.INTER_LINEAR) for c in range(3)], -1)
    vv = cv2.remap(v, mx, my, cv2.INTER_LINEAR)
    out[vv < 0.999] = np.nan
    return out


class Umbilicus:
    def __init__(self, path, scale=1.0):
        j = json.load(open(path))
        pts = j.get('control_points') or j.get('points') or j
        if isinstance(pts, dict):
            pts = pts.get('points', [])
        a = np.array([[p['z'], p['y'], p['x']] for p in pts], np.float64) * scale
        a = a[np.argsort(a[:, 0])]
        self.z, self.y, self.x = a[:, 0], a[:, 1], a[:, 2]

    def at(self, z):
        return float(np.interp(z, self.z, self.x)), float(np.interp(z, self.z, self.y))


@njit(cache=True)
def slice_poly(P, z0):
    """Polyline (W,2) x,y of a winding grid at exact z0 (per column, between the two valid rows bracketing z0)."""
    H, W, _ = P.shape
    out = np.full((W, 2), np.nan, np.float32)
    for j in range(W):
        for i in range(H - 1):
            za = P[i, j, 2]
            zb = P[i + 1, j, 2]
            if np.isnan(za) or np.isnan(zb):
                continue
            if (za - z0) * (zb - z0) <= 0.0:
                t = 0.0 if zb == za else (z0 - za) / (zb - za)
                out[j, 0] = P[i, j, 0] * (1 - t) + P[i + 1, j, 0] * t
                out[j, 1] = P[i, j, 1] * (1 - t) + P[i + 1, j, 1] * t
                break
    return out


@njit(cache=True)
def ray_hits(poly, wid, cx, cy, nray, hit_r, hit_w, hit_n, maxseg):
    """Intersect a polyline with nray rays from (cx, cy); append (r, wid) per ray."""
    two_pi = 2.0 * np.pi
    W = poly.shape[0]
    kmax = hit_r.shape[1]
    for j in range(W - 1):
        x0, y0, x1, y1 = poly[j, 0], poly[j, 1], poly[j + 1, 0], poly[j + 1, 1]
        if np.isnan(x0) or np.isnan(x1):
            continue
        if (x1 - x0) ** 2 + (y1 - y0) ** 2 > maxseg * maxseg:
            continue
        a0 = np.arctan2(y0 - cy, x0 - cx)
        a1 = np.arctan2(y1 - cy, x1 - cx)
        r0 = np.hypot(x0 - cx, y0 - cy)
        r1 = np.hypot(x1 - cx, y1 - cy)
        da = a1 - a0
        if da > np.pi:
            da -= two_pi
        elif da < -np.pi:
            da += two_pi
        if abs(da) < 1e-9:
            continue
        lo = min(a0, a0 + da)
        hi = max(a0, a0 + da)
        m0 = int(np.ceil(lo / two_pi * nray))
        m1 = int(np.floor(hi / two_pi * nray))
        for m in range(m0, m1 + 1):
            th = m * two_pi / nray
            t = (th - a0) / da
            if t < 0.0 or t > 1.0:
                continue
            mm = m % nray
            n = hit_n[mm]
            if n < kmax:
                hit_r[mm, n] = r0 * (1 - t) + r1 * t
                hit_w[mm, n] = wid
                hit_n[mm] = n + 1


@njit(cache=True)
def sheet_candidates(prof, step, thr, minrun, cand, cand_n):
    """Runs of prof > thr (>= minrun samples) -> value-weighted centre radius per ray."""
    nray, T = prof.shape
    cmax = cand.shape[1]
    for m in range(nray):
        n = 0
        t = 0
        while t < T:
            if prof[m, t] > thr:
                s = t
                acc = 0.0
                wsum = 0.0
                while t < T and prof[m, t] > thr:
                    acc += prof[m, t] * t
                    wsum += prof[m, t]
                    t += 1
                if t - s >= minrun and n < cmax:
                    cand[m, n] = (acc / wsum) * step
                    n += 1
            else:
                t += 1
        cand_n[m] = n


@njit(cache=True)
def align_ray(t, rho, sigma_f, dmax_f, cw, cs, naive, out_rho):
    """Order-preserving one-to-one alignment of target radii t (sorted) to sheet radii rho (sorted).
    out_rho[i] = matched sheet radius or NaN."""
    n = t.shape[0]
    m = rho.shape[0]
    for i in range(n):
        out_rho[i] = np.nan
    if n == 0 or m == 0:
        return
    p = np.empty(n)             # local winding spacing: mean of the two neighbour gaps (one at the ends)
    for i in range(n):
        g1 = t[i] - t[i - 1] if i > 0 else -1.0
        g2 = t[i + 1] - t[i] if i < n - 1 else -1.0
        if g1 > 0 and g2 > 0:
            p[i] = 0.5 * (g1 + g2)
        elif g1 > 0:
            p[i] = g1
        elif g2 > 0:
            p[i] = g2
        else:
            p[i] = 20.0
        if p[i] < 4.0:
            p[i] = 4.0
    if naive:
        for i in range(n):
            best = -1
            bd = 1e18
            for j in range(m):
                d = abs(rho[j] - t[i])
                if d < bd:
                    bd = d
                    best = j
            if best >= 0 and bd <= dmax_f * p[i]:
                out_rho[i] = rho[best]
        return
    INF = 1e30
    D = np.full((n + 1, m + 1), INF)
    B = np.zeros((n + 1, m + 1), np.int8)   # 1 match, 2 skip winding, 3 skip sheet
    D[0, 0] = 0.0
    for j in range(1, m + 1):
        D[0, j] = D[0, j - 1] + cs
        B[0, j] = 3
    for i in range(1, n + 1):
        D[i, 0] = D[i - 1, 0] + cw
        B[i, 0] = 2
        for j in range(1, m + 1):
            best = D[i - 1, j] + cw
            b = 2
            v = D[i, j - 1] + cs
            if v < best:
                best = v
                b = 3
            d = abs(t[i - 1] - rho[j - 1])
            if d <= dmax_f * p[i - 1]:
                v = D[i - 1, j - 1] + (d / (sigma_f * p[i - 1])) ** 2
                if v < best:
                    best = v
                    b = 1
            D[i, j] = best
            B[i, j] = b
    i, j = n, m
    while i > 0 or j > 0:
        b = B[i, j]
        if b == 1:
            out_rho[i - 1] = rho[j - 1]
            i -= 1
            j -= 1
        elif b == 2:
            i -= 1
        else:
            j -= 1


@njit(cache=True, parallel=True)
def solve_slice(hit_r, hit_w, hit_o, hit_n, cand, cand_n, sigma_f, dmax_f, cw, cs, naive, res_w, res_d, res_n):
    """Per ray: targets = fitted radius + prior offset, sorted; align; displacement relative to the FITTED radius."""
    nray = hit_r.shape[0]
    for m in prange(nray):
        n = hit_n[m]
        res_n[m] = n
        if n == 0:
            continue
        tgt = hit_r[m, :n] + hit_o[m, :n]
        order = np.argsort(tgt)
        t = tgt[order]
        rho = cand[m, :cand_n[m]]
        out = np.empty(n)
        align_ray(t, rho, sigma_f, dmax_f, cw, cs, naive, out)
        for i in range(n):
            k = order[i]
            res_w[m, i] = hit_w[m, k]
            res_d[m, i] = out[i] - hit_r[m, k]


# --------------------------------------------------------------------------------------------- field smoothing


@njit(cache=True, parallel=True)
def nan_median_filter(F, hk, hm, min_frac):
    """Median over the finite values of a (2hk+1) x (2hm+1) window (slice axis clamped, ray axis wrapped);
    NaN where fewer than min_frac of the window is finite. A median of matched moves is itself a matched move
    (it picks one sheet), unlike a mean, which lands between two sheets where the choice flips."""
    K, M = F.shape
    out = np.full((K, M), np.nan, np.float32)
    need = int(min_frac * (2 * hk + 1) * (2 * hm + 1))
    for k in prange(K):
        buf = np.empty((2 * hk + 1) * (2 * hm + 1), np.float32)
        for m in range(M):
            n = 0
            for dk in range(-hk, hk + 1):
                kk = min(max(k + dk, 0), K - 1)
                for dm in range(-hm, hm + 1):
                    v = F[kk, (m + dm) % M]
                    if not np.isnan(v):
                        buf[n] = v
                        n += 1
            if n >= max(need, 1):
                s = np.sort(buf[:n])
                out[k, m] = s[n // 2]
    return out


def smooth_field(F, med, sig):
    """Displacement prior for one winding over (slices, rays): NaN-aware median of the pass-1 moves where the window
    has support; remaining gaps filled by a normalised Gaussian of that median field (rays wrap around)."""
    v = np.isfinite(F)
    if not v.any():
        return F.copy()
    Mf = nan_median_filter(F.astype(np.float32), med[0] // 2, med[1] // 2, 0.3)
    vm = np.isfinite(Mf)
    sup = gaussian_filter(vm.astype(np.float32), sig, mode=('nearest', 'wrap'))
    num = gaussian_filter(np.where(vm, Mf, 0.0).astype(np.float32), sig, mode=('nearest', 'wrap'))
    fill = np.where(sup > 0.05, num / np.maximum(sup, 1e-6), np.nan)
    return np.where(vm, Mf, fill).astype(np.float32)


# --------------------------------------------------------------------------------------------- main


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--windings', required=True)
    ap.add_argument('--surface', required=True, help='zarr array path (e.g. .../pred.zarr/2)')
    ap.add_argument('--umbilicus', required=True)
    ap.add_argument('--umbilicus-scale', type=float, default=1.0)
    ap.add_argument('--out', required=True)
    ap.add_argument('--z0', type=int, required=True)
    ap.add_argument('--z1', type=int, required=True)
    ap.add_argument('--dz', type=int, default=4)
    ap.add_argument('--rays', type=int, default=2048)
    ap.add_argument('--rmax', type=float, default=4000.0)
    ap.add_argument('--step', type=float, default=0.5, help='radial sampling step (voxels)')
    ap.add_argument('--thr', type=float, default=127.0)
    ap.add_argument('--minrun', type=int, default=1)
    ap.add_argument('--densify', type=int, default=4)
    ap.add_argument('--sigma', type=float, default=0.35, help='match cost scale, fraction of local spacing')
    ap.add_argument('--dmax', type=float, default=0.6, help='pass-1 largest move, fraction of local spacing')
    ap.add_argument('--dmax2', type=float, default=0.45, help='pass-2 largest move from the smoothed target')
    ap.add_argument('--cs', type=float, default=1.0, help='cost of a sheet left without a winding')
    ap.add_argument('--mode', choices=['consistent', 'naive'], default='consistent')
    ap.add_argument('--med', type=int, nargs=2, default=[3, 9], help='median window (slices, rays)')
    ap.add_argument('--sig', type=float, nargs=2, default=[1.5, 6.0], help='Gaussian sigma (slices, rays)')
    ap.add_argument('--bbox', type=int, nargs=4, default=None, help='y0 y1 x0 x1 of the surface data to read')
    ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()
    import numba
    numba.set_num_threads(a.threads)
    t0 = time.time()
    cw = (a.dmax / a.sigma) ** 2           # leaving a winding unmatched costs as much as the largest allowed move
    cw2 = (a.dmax2 / a.sigma) ** 2
    S = zarr.open(a.surface, mode='r')
    umb = Umbilicus(a.umbilicus, a.umbilicus_scale)
    wind = load_windings(a.windings)
    wids = sorted(wind)
    idx = {w: i for i, w in enumerate(wids)}
    print(f'{len(wids)} windings, mode {a.mode}', flush=True)
    zs = np.arange(a.z0, a.z1, a.dz)
    K, M, NW = len(zs), a.rays, len(wids)
    T = int(a.rmax / a.step)
    tt = np.arange(T, dtype=np.float32) * a.step
    th = np.arange(M) * 2 * np.pi / M
    cos, sin = np.cos(th).astype(np.float32), np.sin(th).astype(np.float32)
    y0, y1, x0, x1 = a.bbox if a.bbox else (0, S.shape[1], 0, S.shape[2])
    slab = {'c': None, 'a': None}
    czc = S.chunks[0]

    def get_slice(z):
        c = int(z) // czc
        if slab['c'] != c:          # one chunk row of z at a time: each chunk is decompressed once per pass
            slab['a'] = None
            slab['a'] = np.asarray(S[c * czc:min((c + 1) * czc, S.shape[0]), y0:y1, x0:x1])
            slab['c'] = c
        return slab['a'][int(z) - c * czc].astype(np.float32)

    def run_pass(prior, dmax, cwx, tag):
        F = np.full((NW, K, M), np.nan, np.float32)
        P = np.zeros((NW, K, M), bool)
        st = dict(rays=0, hits=0, matched=0, sheets=0, fold_hits=0)
        slab['c'] = None
        for k, z in enumerate(zs):
            sl = get_slice(z)
            cx, cy = umb.at(z)
            mx = (cx - x0 + tt[None, :] * cos[:, None]).astype(np.float32)
            my = (cy - y0 + tt[None, :] * sin[:, None]).astype(np.float32)
            prof = cv2.remap(sl, mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT, borderValue=0)
            cand = np.zeros((M, 1024), np.float64)
            cand_n = np.zeros(M, np.int64)
            sheet_candidates(prof, a.step, a.thr, a.minrun, cand, cand_n)
            kmax = 2 * NW + 16
            hit_r = np.zeros((M, kmax), np.float64)
            hit_w = np.zeros((M, kmax), np.int64)
            hit_n = np.zeros(M, np.int64)
            for w in wids:
                poly = slice_poly(wind[w][1], float(z))
                ray_hits(poly, w, cx, cy, M, hit_r, hit_w, hit_n, 120.0)
            hit_o = np.zeros((M, kmax), np.float64)
            if prior is not None:
                mm = np.repeat(np.arange(M)[:, None], kmax, 1)
                wi = np.vectorize(lambda q: idx.get(int(q), 0))(hit_w)
                pv = prior[wi, k, mm]
                hit_o = np.where(np.isfinite(pv), pv, 0.0).astype(np.float64)
            res_w = np.zeros((M, kmax), np.int64)
            res_d = np.full((M, kmax), np.nan)
            res_n = np.zeros(M, np.int64)
            solve_slice(hit_r, hit_w, hit_o, hit_n, cand, cand_n, a.sigma, dmax, cwx, a.cs, a.mode == 'naive',
                        res_w, res_d, res_n)
            for m in range(M):
                n = res_n[m]
                for i in range(n):
                    wi = idx[int(res_w[m, i])]
                    if P[wi, k, m]:           # a winding crossing one ray twice (fold): keep neither
                        F[wi, k, m] = np.nan
                        st['fold_hits'] += 1
                        continue
                    P[wi, k, m] = True
                    F[wi, k, m] = res_d[m, i]
                st['rays'] += 1
                st['hits'] += int(n)
                st['matched'] += int(np.isfinite(res_d[m, :n]).sum())
                st['sheets'] += int(cand_n[m])
            if k % 50 == 0:
                print(f'{tag} slice {k + 1}/{K} z {z}: {time.time() - t0:.0f} s', flush=True)
        return F, P, st

    F1, Pres, st1 = run_pass(None, a.dmax, cw, 'pass 1')
    Fs = np.stack([smooth_field(F1[i], tuple(a.med), tuple(a.sig)) for i in range(NW)])
    F2, _, st2 = run_pass(Fs, a.dmax2, cw2, 'pass 2')
    print('pass 1', st1, '\npass 2', st2, flush=True)
    # ---- final displacement per (winding, slice, ray): pass-2 match, else smoothed field
    Ffin = np.where(np.isfinite(F2), F2, Fs)
    Fsrc = np.where(np.isfinite(F2), 0, np.where(np.isfinite(Fs), 1, 3)).astype(np.uint8)
    # ---- move nodes
    os.makedirs(a.out, exist_ok=True)
    rep = dict(args=vars(a), pass1=st1, pass2=st2, windings={})
    cxy_z = np.array([umb.at(float(z)) for z in zs])
    for w in wids:
        name, P0, meta = wind[w]
        P = densify(P0, a.densify)
        wi = idx[w]
        H, W, _ = P.shape
        flag = np.full((H, W), 3, np.uint8)
        dout = np.zeros((H, W), np.float32)
        ok = np.isfinite(P).all(-1)
        Z = P[..., 2]
        inz = ok & (Z >= zs[0] - a.dz / 2) & (Z <= zs[-1] + a.dz / 2)
        rr, cc = np.nonzero(inz)
        z = Z[rr, cc]
        k = np.clip(np.rint((z - zs[0]) / a.dz).astype(int), 0, K - 1)
        cxy = cxy_z[k]
        dx, dy = P[rr, cc, 0] - cxy[:, 0], P[rr, cc, 1] - cxy[:, 1]
        ang = np.mod(np.arctan2(dy, dx), 2 * np.pi)
        m = np.rint(ang / (2 * np.pi) * M).astype(int) % M
        d = Ffin[wi, k, m]
        src = Fsrc[wi, k, m]
        good = np.isfinite(d)
        flag[rr, cc] = np.where(good, src, 3)
        r = np.hypot(dx, dy)
        scale = np.where(r > 1e-6, (r + np.nan_to_num(d)) / np.maximum(r, 1e-6), 1.0)
        Pn = P.copy()
        Pn[rr[good], cc[good], 0] = (cxy[:, 0] + dx * scale)[good]
        Pn[rr[good], cc[good], 1] = (cxy[:, 1] + dy * scale)[good]
        dout[rr[good], cc[good]] = d[good]
        od = os.path.join(a.out, name)
        os.makedirs(od, exist_ok=True)
        for i, c in enumerate('xyz'):
            tifffile.imwrite(os.path.join(od, f'{c}.tif'), np.where(ok, Pn[..., i], -1).astype(np.float32))
        tifffile.imwrite(os.path.join(od, 'snap_flag.tif'), flag)
        tifffile.imwrite(os.path.join(od, 'snap_d.tif'), dout)
        meta = dict(meta)
        sc = meta.get('scale', [0.05, 0.05])
        meta['scale'] = [sc[0] * a.densify, sc[1] * a.densify]
        meta['source'] = f"consistent_snap ({a.mode}) of {meta.get('source', 'fit')}"
        meta['consistent_snap'] = dict(mode=a.mode, surface=a.surface, sigma=a.sigma, dmax=a.dmax, dmax2=a.dmax2,
                                       cs=a.cs, densify=a.densify)
        json.dump(meta, open(os.path.join(od, 'meta.json'), 'w'), indent=1)
        fl = flag[rr, cc]
        rep['windings'][name] = dict(nodes=int(ok.sum()), in_z=int(inz.sum()), matched=int((fl == 0).sum()),
                                     smoothed=int((fl == 1).sum()), unmoved=int((fl == 3).sum()),
                                     median_abs_move=float(np.median(np.abs(d[good]))) if good.any() else None)
    tot = {k2: sum(v[k2] for v in rep['windings'].values()) for k2 in ('nodes', 'in_z', 'matched', 'smoothed', 'unmoved')}
    rep['total'] = tot
    rep['seconds'] = round(time.time() - t0, 1)
    json.dump(rep, open(os.path.join(a.out, 'report.json'), 'w'), indent=1)
    print('DONE', json.dumps(tot), f'{time.time() - t0:.0f} s', flush=True)


if __name__ == '__main__':
    main()
