"""track_snap.py: lock the windings of a spiral fit onto physical sheets of a surface prediction CONTINUOUSLY, so the
snap does not turn the fit's slow drift across sheets into stair-steps (abrupt sheet jumps).

Why: a fit whose local slope differs from the sheets crosses them diagonally. A per-node snap (nearest sheet, or the
order-consistent per-ray alignment of consistent_snap.py) then jumps from sheet j to sheet j+1 wherever the fit has
drifted half a pitch: stair-steps, i.e. sheet switches in cross-section.

How: one joint sweep over the band. Slices are visited outward from the middle slice (up, then down); in each slice the
rays are visited in angular order from ray m0. On every ray all windings present are assigned to predicted sheets at
once by the exact order-preserving one-to-one alignment of consistent_snap.py (no two windings on one sheet, no
crossings), but the target of each winding is where its sheet CONTINUES:
    z-carry:   its displacement (sheet radius - fitted radius) on the same ray in the previous visited slice, else
    ray-carry: its displacement on the previous ray of this slice (linked by fitted radius, so a track crosses the
               fit's winding seam onto the neighbouring winding id), else
    acquire:   its fitted radius (wider tolerance; counted).
Tracked targets use a tight tolerance (--sig-t, --dmax-t, fractions of the local pitch), so a winding keeps its sheet and
its displacement may drift away from the fit by up to --reach pitches. A carried displacement survives --carry visited
cells without a match (gaps in the prediction), then the winding re-acquires.

Nodes (grid densified f x) move radially: matched cells (flag 0), gaps filled from a robust NaN-median field (flag 1),
else unmoved (flag 3). Output as consistent_snap.py: tifxyz + snap_flag.tif + snap_d.tif per winding, report.json.
CPU only.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import cv2
import numpy as np
import tifffile
import zarr
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consistent_snap import (Umbilicus, densify, load_windings, ray_hits, sheet_candidates, slice_poly,  # noqa: E402
                             smooth_field)


@njit(cache=True)
def fitted_table(hit_r, hit_w, hit_n, widx, T, Pch):
    """Per (winding, ray): fitted radius T (NaN if none or folded) and local pitch Pch from the radial neighbours."""
    M = hit_r.shape[0]
    NW = T.shape[0]
    fold = np.zeros((NW, M), np.bool_)
    for m in range(M):
        n = hit_n[m]
        if n == 0:
            continue
        order = np.argsort(hit_r[m, :n])
        rs = hit_r[m, :n][order]
        ws = hit_w[m, :n][order]
        for i in range(n):
            wi = widx[ws[i]]
            if wi < 0:
                continue
            if not np.isnan(T[wi, m]) or fold[wi, m]:
                fold[wi, m] = True           # second crossing of one ray: a fold, keep neither
                T[wi, m] = np.nan
                continue
            T[wi, m] = rs[i]
            lo = rs[i] - rs[i - 1] if i > 0 else np.nan
            hi = rs[i + 1] - rs[i] if i + 1 < n else np.nan
            if np.isnan(lo):
                p = hi
            elif np.isnan(hi):
                p = lo
            else:
                p = 0.5 * (lo + hi)
            Pch[wi, m] = p


@njit(cache=True)
def align_ray2(t, p, rho, sig, dmax, cs, out_rho):
    """align_ray of consistent_snap.py with per-winding pitch p, tolerance sig and largest move dmax (fractions of p);
    leaving winding i unmatched costs (dmax_i / sig_i)^2, a sheet without a winding costs cs."""
    n = t.shape[0]
    m = rho.shape[0]
    for i in range(n):
        out_rho[i] = np.nan
    if n == 0 or m == 0:
        return
    INF = 1e30
    D = np.full((n + 1, m + 1), INF)
    B = np.zeros((n + 1, m + 1), np.int8)
    D[0, 0] = 0.0
    for j in range(1, m + 1):
        D[0, j] = D[0, j - 1] + cs
        B[0, j] = 3
    for i in range(1, n + 1):
        cw = (dmax[i - 1] / sig[i - 1]) ** 2
        D[i, 0] = D[i - 1, 0] + cw
        B[i, 0] = 2
        lim = dmax[i - 1] * p[i - 1]
        sp = sig[i - 1] * p[i - 1]
        for j in range(1, m + 1):
            best = D[i - 1, j] + cw
            b = 2
            v = D[i, j - 1] + cs
            if v < best:
                best = v
                b = 3
            d = abs(t[i - 1] - rho[j - 1])
            if d <= lim:
                v = D[i - 1, j - 1] + (d / sp) ** 2
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


@njit(cache=True)
def sweep(T, Pch, cand, cand_n, ks, m0, sig_a, dmax_a, sig_t, dmax_t, cs, carry, reach, Zc, Zm, Zt, D, FL, stats):
    """Visit slices ks in order (z-carry state Zc/Zm/Zt (NW, M) carried between them) and rays from m0. Writes D (K, NW, M)
    displacement or NaN and FL (K, NW, M) mode: 0 ray-carry, 1 z-carry, 2 acquire, 255 none.
    stats: [cells, matched, acquired, ray_tracked, z_tracked, over_reach]"""
    K, NW, M = T.shape
    for kk in range(len(ks)):
        k = ks[kk]
        car = np.full(NW, np.nan)
        mis = np.zeros(NW, np.int64)
        prv = np.full(NW, np.nan)
        for step in range(M):
            m = (m0 + step) % M
            ws0 = np.nonzero(~np.isnan(T[k, :, m]))[0]
            n = len(ws0)
            if n == 0:
                continue
            ws = ws0[np.argsort(T[k, ws0, m])]
            q = np.empty(n)
            sg = np.empty(n)
            dmx = np.empty(n)
            pp = np.empty(n)
            mode = np.zeros(n, np.int64)
            cmis = np.zeros(n, np.int64)
            for i in range(n):
                w = ws[i]
                t = T[k, w, m]
                p = Pch[k, w, m]
                if np.isnan(p):
                    p = 12.0
                if p < 4.0:
                    p = 4.0
                pp[i] = p
                c = np.nan
                md = 2
                if not np.isnan(Zc[w, m]) and Zm[w, m] < carry and abs(Zt[w, m] - t) < 0.5 * p:
                    c = Zc[w, m]
                    md = 1
                    cmis[i] = Zm[w, m]
                else:
                    src = -1
                    if not np.isnan(prv[w]) and abs(prv[w] - t) < 0.5 * p:
                        src = w
                    else:
                        bd = 0.5 * p
                        for v in range(NW):
                            if not np.isnan(prv[v]):
                                dd = abs(prv[v] - t)
                                if dd < bd:
                                    bd = dd
                                    src = v
                    if src >= 0 and not np.isnan(car[src]) and mis[src] < carry:
                        c = car[src]
                        md = 0
                        cmis[i] = mis[src]
                if md == 2:
                    q[i] = t
                    sg[i] = sig_a
                    dmx[i] = dmax_a
                else:
                    q[i] = t + c
                    sg[i] = sig_t
                    dmx[i] = dmax_t
                mode[i] = md
            order = np.argsort(q)
            out = np.empty(n)
            rho = cand[k, m, :cand_n[k, m]].astype(np.float64)
            align_ray2(q[order], pp[order], rho, sg[order], dmx[order], cs, out)
            for ii in range(n):
                i = order[ii]
                w = ws[i]
                t = T[k, w, m]
                r = out[ii]
                stats[0] += 1
                prv[w] = t
                if not np.isnan(r) and abs(r - t) <= reach * pp[i]:
                    d = r - t
                    D[k, w, m] = d
                    FL[k, w, m] = mode[i]
                    car[w] = d
                    mis[w] = 0
                    Zc[w, m] = d
                    Zm[w, m] = 0
                    Zt[w, m] = t
                    stats[1] += 1
                    if mode[i] == 2:
                        stats[2] += 1
                    elif mode[i] == 0:
                        stats[3] += 1
                    else:
                        stats[4] += 1
                else:
                    if not np.isnan(r):
                        stats[5] += 1
                    if mode[i] == 2:
                        car[w] = np.nan
                        mis[w] = 0
                    else:
                        car[w] = q[i] - t
                        mis[w] = cmis[i] + 1
                    Zm[w, m] += 1
                    Zt[w, m] = t


@njit(cache=True)
def _hpush(hp, hk, hm, n, p, k, m):
    i = n
    hp[i] = p
    hk[i] = k
    hm[i] = m
    while i > 0:
        par = (i - 1) // 2
        if hp[par] >= hp[i]:
            break
        hp[par], hp[i] = hp[i], hp[par]
        hk[par], hk[i] = hk[i], hk[par]
        hm[par], hm[i] = hm[i], hm[par]
        i = par
    return n + 1


@njit(cache=True)
def _hpop(hp, hk, hm, n):
    p, k, m = hp[0], hk[0], hm[0]
    n -= 1
    hp[0], hk[0], hm[0] = hp[n], hk[n], hm[n]
    i = 0
    while True:
        lft = 2 * i + 1
        rgt = lft + 1
        big = i
        if lft < n and hp[lft] > hp[big]:
            big = lft
        if rgt < n and hp[rgt] > hp[big]:
            big = rgt
        if big == i:
            break
        hp[big], hp[i] = hp[i], hp[big]
        hk[big], hk[i] = hk[i], hk[big]
        hm[big], hm[i] = hm[i], hm[big]
        i = big
    return p, k, m, n


@njit(cache=True)
def _link(T, k2, m2, w, t, p):
    """Winding at ray (k2, m2) that continues winding w (fitted radius t): the same id if its fitted radius is within half
    a pitch, else the nearest fitted radius within half a pitch (crosses the fit's winding seam), else -1."""
    v = T[k2, w, m2]
    if not np.isnan(v) and abs(v - t) < 0.5 * p:
        return w
    best = 0.5 * p
    bv = -1
    for u in range(T.shape[1]):
        x = T[k2, u, m2]
        if not np.isnan(x):
            dd = abs(x - t)
            if dd < best:
                best = dd
                bv = u
    return bv


@njit(cache=True)
def _solve_ray(T, Pch, cand, cand_n, k, m, visited, Q, Dd, Cc, Mi, FL, sig_a, dmax_a, sig_t, dmax_t, cs, carry, reach,
               stats):
    K, NW, M = T.shape
    ws0 = np.nonzero(~np.isnan(T[k, :, m]))[0]
    n = len(ws0)
    if n == 0:
        Q[k, m] = 0.0
        return
    ws = ws0[np.argsort(T[k, ws0, m])]
    q = np.empty(n)
    sg = np.empty(n)
    dmx = np.empty(n)
    pp = np.empty(n)
    mode = np.zeros(n, np.int64)
    cmis = np.zeros(n, np.int64)
    nk = np.array([k - 1, k + 1, k, k])
    nm = np.array([m, m, (m - 1) % M, (m + 1) % M])
    vals = np.empty(4)
    qual = np.empty(4)
    miss = np.empty(4, np.int64)
    for i in range(n):
        w = ws[i]
        t = T[k, w, m]
        p = Pch[k, w, m]
        if np.isnan(p):
            p = 12.0
        if p < 4.0:
            p = 4.0
        pp[i] = p
        nv = 0
        for j in range(4):
            k2 = nk[j]
            m2 = nm[j]
            if k2 < 0 or k2 >= K:
                continue
            if not visited[k2, m2]:
                continue
            v = _link(T, k2, m2, w, t, p)
            if v < 0:
                continue
            c = Cc[k2, v, m2]
            if np.isnan(c) or Mi[k2, v, m2] >= carry:
                continue
            vals[nv] = c
            qual[nv] = Q[k2, m2]
            miss[nv] = Mi[k2, v, m2]
            nv += 1
        if nv == 0:
            q[i] = t
            sg[i] = sig_a
            dmx[i] = dmax_a
            mode[i] = 2
            continue
        vmin = vals[0]
        vmax = vals[0]
        for j in range(1, nv):
            vmin = min(vmin, vals[j])
            vmax = max(vmax, vals[j])
        if vmax - vmin <= 0.5 * p:
            s = np.sort(vals[:nv])
            if nv % 2 == 1:
                c = s[nv // 2]
            else:
                c = 0.5 * (s[nv // 2 - 1] + s[nv // 2])
            mm = miss[0]
            for j in range(1, nv):
                mm = min(mm, miss[j])
        else:
            jb = 0
            for j in range(1, nv):
                if qual[j] > qual[jb]:
                    jb = j
            c = vals[jb]
            mm = miss[jb]
            stats[6] += 1
        q[i] = t + c
        sg[i] = sig_t
        dmx[i] = dmax_t
        mode[i] = 0
        cmis[i] = mm
    order = np.argsort(q)
    out = np.empty(n)
    rho = cand[k, m, :cand_n[k, m]].astype(np.float64)
    align_ray2(q[order], pp[order], rho, sg[order], dmx[order], cs, out)
    good = 0
    for ii in range(n):
        i = order[ii]
        w = ws[i]
        t = T[k, w, m]
        r = out[ii]
        stats[0] += 1
        if not np.isnan(r) and abs(r - t) <= reach * pp[i]:
            Dd[k, w, m] = r - t
            Cc[k, w, m] = r - t
            Mi[k, w, m] = 0
            FL[k, w, m] = mode[i]
            good += 1
            stats[1] += 1
            if mode[i] == 2:
                stats[2] += 1
        else:
            if not np.isnan(r):
                stats[5] += 1
            if mode[i] == 2:
                Cc[k, w, m] = np.nan
                Mi[k, w, m] = 0
            else:
                Cc[k, w, m] = q[i] - t
                Mi[k, w, m] = cmis[i] + 1
    Q[k, m] = good / n


@njit(cache=True)
def grow(T, Pch, cand, cand_n, k0, m0, sig_a, dmax_a, sig_t, dmax_t, cs, carry, reach, Dd, Cc, Mi, FL, Q, order_k,
         order_m, stats):
    """Quality-guided growth over the (slice, ray) cylinder from the seed ray (k0, m0): the unsolved neighbour of the best
    solved ray is solved next, its targets carried from all solved neighbours (median; best-quality neighbour where they
    disagree by more than half a pitch, counted in stats[6]). stats: [cells, matched, acquired, -, -, over_reach,
    disagreements]."""
    K, NW, M = T.shape
    visited = np.zeros((K, M), np.bool_)
    cap = 4 * K * M + 16
    hp = np.empty(cap)
    hk = np.empty(cap, np.int64)
    hm = np.empty(cap, np.int64)
    n = _hpush(hp, hk, hm, 0, 1.0, k0, m0)
    cnt = 0
    while n > 0:
        p, k, m, n = _hpop(hp, hk, hm, n)
        if visited[k, m]:
            continue
        _solve_ray(T, Pch, cand, cand_n, k, m, visited, Q, Dd, Cc, Mi, FL, sig_a, dmax_a, sig_t, dmax_t, cs, carry,
                   reach, stats)
        visited[k, m] = True
        order_k[cnt] = k
        order_m[cnt] = m
        cnt += 1
        qq = Q[k, m]
        for j in range(4):
            if j == 0:
                k2 = k - 1
                m2 = m
            elif j == 1:
                k2 = k + 1
                m2 = m
            elif j == 2:
                k2 = k
                m2 = (m - 1) % M
            else:
                k2 = k
                m2 = (m + 1) % M
            if k2 >= 0 and k2 < K:
                if not visited[k2, m2]:
                    n = _hpush(hp, hk, hm, n, qq, k2, m2)
    return cnt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--windings', required=True)
    ap.add_argument('--surface', required=True)
    ap.add_argument('--umbilicus', required=True)
    ap.add_argument('--umbilicus-scale', type=float, default=1.0)
    ap.add_argument('--out', required=True)
    ap.add_argument('--z0', type=int, required=True)
    ap.add_argument('--z1', type=int, required=True)
    ap.add_argument('--dz', type=int, default=4)
    ap.add_argument('--rays', type=int, default=2048)
    ap.add_argument('--rmax', type=float, default=4000.0)
    ap.add_argument('--step', type=float, default=0.5)
    ap.add_argument('--thr', type=float, default=127.0)
    ap.add_argument('--minrun', type=int, default=1)
    ap.add_argument('--densify', type=int, default=2)
    ap.add_argument('--sig-a', type=float, default=0.35, help='acquire tolerance, fraction of local pitch')
    ap.add_argument('--dmax-a', type=float, default=0.6, help='acquire largest move, fraction of local pitch')
    ap.add_argument('--sig-t', type=float, default=0.15, help='tracking tolerance, fraction of local pitch')
    ap.add_argument('--dmax-t', type=float, default=0.35, help='tracking largest step, fraction of local pitch')
    ap.add_argument('--cs', type=float, default=1.0, help='cost of a sheet left without a winding')
    ap.add_argument('--carry', type=int, default=16, help='cells a carried displacement survives without a match')
    ap.add_argument('--reach', type=float, default=2.5, help='largest total move from the fit, in local pitches')
    ap.add_argument('--m0', type=int, default=0, help='first ray of every slice sweep (sweep mode)')
    ap.add_argument('--mode', choices=['grow', 'sweep'], default='grow')
    ap.add_argument('--med', type=int, nargs=2, default=[3, 9])
    ap.add_argument('--sig', type=float, nargs=2, default=[1.5, 6.0])
    ap.add_argument('--bbox', type=int, nargs=4, default=None)
    ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()
    import numba
    numba.set_num_threads(a.threads)
    t0 = time.time()
    S = zarr.open(a.surface, mode='r')
    umb = Umbilicus(a.umbilicus, a.umbilicus_scale)
    wind = load_windings(a.windings)
    wids = sorted(wind)
    widx = np.full(max(wids) + 1, -1, np.int64)
    for i, w in enumerate(wids):
        widx[w] = i
    print(f'{len(wids)} windings, {a.mode} mode', flush=True)
    zs = np.arange(a.z0, a.z1, a.dz)
    K, M, NW = len(zs), a.rays, len(wids)
    Tn = int(a.rmax / a.step)
    tt = np.arange(Tn, dtype=np.float32) * a.step
    th = np.arange(M) * 2 * np.pi / M
    cos, sin = np.cos(th).astype(np.float32), np.sin(th).astype(np.float32)
    y0, y1, x0, x1 = a.bbox if a.bbox else (0, S.shape[1], 0, S.shape[2])
    slab = {'c': None, 'a': None}
    czc = S.chunks[0]

    def get_slice(z):
        c = int(z) // czc
        if slab['c'] != c:
            slab['a'] = None
            slab['a'] = np.asarray(S[c * czc:min((c + 1) * czc, S.shape[0]), y0:y1, x0:x1])
            slab['c'] = c
        return slab['a'][int(z) - c * czc].astype(np.float32)

    T = np.full((K, NW, M), np.nan)
    Pch = np.full((K, NW, M), np.nan)
    cands, cns = [], []
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
        fitted_table(hit_r, hit_w, hit_n, widx, T[k], Pch[k])
        cands.append(cand[:, :max(int(cand_n.max()), 1)].astype(np.float32))
        cns.append(cand_n)
        if k % 50 == 0:
            print(f'prep slice {k + 1}/{K} z {z}: {time.time() - t0:.0f} s', flush=True)
    C = max(c.shape[1] for c in cands)
    CAND = np.zeros((K, M, C), np.float32)
    for k, c in enumerate(cands):
        CAND[k, :, :c.shape[1]] = c
    CN = np.stack(cns)
    del cands, cns
    D = np.full((K, NW, M), np.nan)
    FL = np.full((K, NW, M), 255, np.uint8)
    kmid = K // 2
    if a.mode == 'sweep':
        stats = np.zeros(6, np.int64)
        Zc = np.full((NW, M), np.nan)
        Zm = np.zeros((NW, M), np.int64)
        Zt = np.full((NW, M), np.nan)
        sweep(T, Pch, CAND, CN, np.array([kmid], np.int64), a.m0, a.sig_a, a.dmax_a, a.sig_t, a.dmax_t, a.cs,
              a.carry, a.reach, Zc, Zm, Zt, D, FL, stats)
        Zc0, Zm0, Zt0 = Zc.copy(), Zm.copy(), Zt.copy()
        sweep(T, Pch, CAND, CN, np.arange(kmid + 1, K, dtype=np.int64), a.m0, a.sig_a, a.dmax_a, a.sig_t, a.dmax_t,
              a.cs, a.carry, a.reach, Zc, Zm, Zt, D, FL, stats)
        sweep(T, Pch, CAND, CN, np.arange(kmid - 1, -1, -1, dtype=np.int64), a.m0, a.sig_a, a.dmax_a, a.sig_t,
              a.dmax_t, a.cs, a.carry, a.reach, Zc0, Zm0, Zt0, D, FL, stats)
        st = dict(zip(['cells', 'matched', 'acquired', 'ray_tracked', 'z_tracked', 'over_reach'], map(int, stats)))
        Fcar = D
    else:
        # seed: the ray of the middle slice where an acquire-only alignment matches the largest share of windings
        stats = np.zeros(7, np.int64)
        Cc = np.full((K, NW, M), np.nan)
        Mi = np.zeros((K, NW, M), np.int16)
        Q = np.zeros((K, M))
        vis0 = np.zeros((K, M), np.bool_)
        s0 = np.zeros(7, np.int64)
        Dt, Ct, Mt, Ft = D.copy(), Cc.copy(), Mi.copy(), FL.copy()
        for m in range(0, M, 8):
            _solve_ray(T, Pch, CAND, CN, kmid, m, vis0, Q, Dt, Ct, Mt, Ft, a.sig_a, a.dmax_a, a.sig_t, a.dmax_t, a.cs,
                       a.carry, a.reach, s0)
        mseed = int(np.argmax(Q[kmid]))
        print(f'seed slice {kmid} ray {mseed}: acquire share {Q[kmid, mseed]:.2f}', flush=True)
        del Dt, Ct, Mt, Ft
        Q[:] = 0.0
        ok_ = np.zeros(K * M, np.int64)
        om_ = np.zeros(K * M, np.int64)
        cnt = grow(T, Pch, CAND, CN, kmid, mseed, a.sig_a, a.dmax_a, a.sig_t, a.dmax_t, a.cs, a.carry, a.reach, D, Cc,
                   Mi, FL, Q, ok_, om_, stats)
        st = dict(zip(['cells', 'matched', 'acquired', '-', '-', 'over_reach', 'disagreements'], map(int, stats)))
        st['rays_solved'] = int(cnt)
        st['seed'] = [int(kmid), mseed]
        Fcar = np.where(np.isfinite(D), D, np.where(Mi < a.carry, Cc, np.nan))
    print(a.mode, st, f'{time.time() - t0:.0f} s', flush=True)
    # per winding field over (slice, ray) for the node moves; gaps filled from the robust median field
    F = np.transpose(D, (1, 0, 2)).astype(np.float32)
    if a.mode == 'sweep':
        Fs = np.stack([smooth_field(F[i], tuple(a.med), tuple(a.sig)) for i in range(NW)])
    else:
        Fs = np.transpose(Fcar, (1, 0, 2)).astype(np.float32)     # carried displacement: continuous with neighbours
    Ffin = np.where(np.isfinite(F), F, Fs)
    Fsrc = np.where(np.isfinite(F), 0, np.where(np.isfinite(Fs), 1, 3)).astype(np.uint8)
    os.makedirs(a.out, exist_ok=True)
    rep = dict(args=vars(a), run=st, windings={})
    cxy_z = np.array([umb.at(float(z)) for z in zs])
    for w in wids:
        name, P0, meta = wind[w]
        P = densify(P0, a.densify)
        wi = widx[w]
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
        meta['source'] = f"track_snap of {meta.get('source', 'fit')}"
        meta['track_snap'] = {k2: v for k2, v in vars(a).items() if k2 not in ('windings', 'out')}
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
