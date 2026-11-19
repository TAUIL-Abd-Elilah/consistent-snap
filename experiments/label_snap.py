"""label_snap.py: snap all windings of a spiral fit onto a surface prediction with ONE winding label per predicted sheet,
so a snapped winding follows its sheet instead of stepping to the next sheet where the fit drifts.

1. Sheet nodes: runs of the prediction (> --thr) on --rays rays from the umbilicus, every --dz slices (as v1).
   phi(node) = the fit's fractional winding id, interpolated between the fitted crossings on that ray.
2. Links: between neighbouring rays (same slice) and neighbouring slices (same ray), mutual nearest sheet nodes within
   --link-dmax x local pitch (1-D mutual nearest pairs are one-to-one and keep the radial order).
   cost = (dr / (--link-sigma x pitch))^2.  Label jump across a link: delta = round(phi(v) - phi(u)) (0 except at the
   fit's winding seam); links with |phi(v) - phi(u) - delta| > --phi-tol or a NaN phi are not used.
3. Union-find with integer potentials (label(x) - label(root)), cheapest links first. A link is refused when, among
   the --win nearest sheet nodes on either of its two rays, the merged labels would not increase strictly outward.
4. Each component gets one integer offset to the fit: round(median(phi - potential)); components smaller than
   --min-comp nodes or without phi are dropped.
5. Per ray, the longest strictly increasing run of labels (radial order) is kept; other nodes are dropped.
6. Fitted winding w on a ray moves radially to the sheet node labelled w on that ray, if within --reach x pitch.
   Every other node is a hole (-1), as in the matched-only arm of consistent_snap.
Output: tifxyz per winding (+ snap_flag.tif 0 = moved / 3 = hole, snap_d.tif) and report.json.  CPU only.
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
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root
from consistent_snap import Umbilicus, densify, load_windings, ray_hits, sheet_candidates, slice_poly  # noqa: E402


# --------------------------------------------------------------------------------------------- per slice: nodes, phi


@njit(cache=True)
def ray_phi(cand, cand_n, hit_r, hit_w, hit_n, R, PHI, PIT, HR, HW, HN):
    """Per ray: drop windings that cross the ray twice (folds), sort crossings, local pitch (median gap, >= 4,
    default 12), phi of each sheet node by linear interpolation of winding id over radius (NaN outside)."""
    M = cand.shape[0]
    for m in range(M):
        n = hit_n[m]
        r = hit_r[m, :n].copy()
        w = hit_w[m, :n].copy()
        keep = np.ones(n, np.bool_)
        for i in range(n):
            for j in range(i + 1, n):
                if w[i] == w[j]:
                    keep[i] = False
                    keep[j] = False
        r = r[keep]
        w = w[keep]
        o = np.argsort(r)
        r = r[o]
        w = w[o]
        nn = len(r)
        HN[m] = nn
        for i in range(nn):
            HR[m, i] = r[i]
            HW[m, i] = w[i]
        if nn >= 2:
            g = np.sort(r[1:] - r[:-1])
            p = g[len(g) // 2]
            PIT[m] = max(p, 4.0)
        else:
            PIT[m] = 12.0
        c = cand_n[m]
        j = 0
        for i in range(c):
            rho = cand[m, i]
            R[m, i] = rho
            PHI[m, i] = np.nan
            if nn < 2:
                continue
            while j < nn - 2 and r[j + 1] <= rho:
                j += 1
            if rho >= r[0] and rho <= r[nn - 1] and r[j + 1] > r[j]:
                t = (rho - r[j]) / (r[j + 1] - r[j])
                PHI[m, i] = w[j] + t * (w[j + 1] - w[j])


# --------------------------------------------------------------------------------------------- links


@njit(cache=True)
def _mutual(ra, rb, tol, sig, ia, ib, out_v, out_c, base_a, base_b):
    """Mutual nearest pairs between sorted radius lists ra, rb within tol. Writes, for node base_a+i, the partner
    base_b+j and cost ((ra-rb)/sig)^2 into out_v/out_c at index ia+i."""
    na = len(ra)
    nb = len(rb)
    if na == 0 or nb == 0:
        return
    nba = np.empty(na, np.int64)
    j = 0
    for i in range(na):
        while j < nb - 1 and abs(rb[j + 1] - ra[i]) <= abs(rb[j] - ra[i]):
            j += 1
        nba[i] = j
    nab = np.empty(nb, np.int64)
    i = 0
    for j in range(nb):
        while i < na - 1 and abs(ra[i + 1] - rb[j]) <= abs(ra[i] - rb[j]):
            i += 1
        nab[j] = i
    for i in range(na):
        j = nba[i]
        if nab[j] == i:
            d = abs(ra[i] - rb[j])
            if d <= tol:
                out_v[ia + i] = base_b + j
                out_c[ia + i] = (d / sig) ** 2


@njit(cache=True, parallel=True)
def build_links(R, S, PIT, K, M, dmax_f, sig_f, E1, C1, E2, C2):
    """E1: ray link (k, m) -> (k, m+1 mod M); E2: slice link (k, m) -> (k+1, m)."""
    for k in prange(K):
        for m in range(M):
            a0 = S[k * M + m]
            a1 = S[k * M + m + 1]
            if a1 == a0:
                continue
            p = PIT[k * M + m]
            ra = R[a0:a1]
            m2 = (m + 1) % M
            b0 = S[k * M + m2]
            b1 = S[k * M + m2 + 1]
            _mutual(ra, R[b0:b1], dmax_f * p, sig_f * p, a0, 0, E1, C1, 0, b0)
            if k + 1 < K:
                b0 = S[(k + 1) * M + m]
                b1 = S[(k + 1) * M + m + 1]
                _mutual(ra, R[b0:b1], dmax_f * p, sig_f * p, a0, 0, E2, C2, 0, b0)


# --------------------------------------------------------------------------------------------- union-find with potentials


@njit(cache=True)
def _find(x, parent, pot):
    r = x
    acc = 0
    while parent[r] != r:
        acc += pot[r]
        r = parent[r]
    y = x
    a = acc
    while parent[y] != y:
        nxt = parent[y]
        py = pot[y]
        parent[y] = r
        pot[y] = a
        a -= py
        y = nxt
    return r, acc


@njit(cache=True)
def _order_ok(ray0, ray1, x0, R, parent, pot, ra, rb, lab_b_shift, win):
    """Nodes of the ray [ray0, ray1) within win positions of x0: those in component ra (label rel ra) or rb (label rel
    rb + lab_b_shift = rel ra) must have labels strictly increasing with radius (nodes are sorted by radius)."""
    lo = max(ray0, x0 - win)
    hi = min(ray1, x0 + win + 1)
    last = -2147483647
    for x in range(lo, hi):
        rx, lx = _find(x, parent, pot)
        if rx == ra:
            lab = lx
        elif rx == rb:
            lab = lx + lab_b_shift
        else:
            continue
        if lab <= last:
            return False
        last = lab
    return True


@njit(cache=True)
def union_all(order, nE1, E1, E2, PHI, R, RAY, S, parent, pot, size, dsum, dcnt, phi_tol, comp_tol, win, stats):
    """stats: [used, merged, cycle_ok, cycle_bad, refused_order, skipped_phi, refused_dev]. dsum/dcnt (per root): sum
    and count of dev = phi - label (rel root); a merge is refused when the two components' mean dev differ by more than
    comp_tol after the label shift, so a component cannot creep across sheets through many small phi changes."""
    for e in order:
        if e < nE1:
            u = e
            v = E1[e]
        else:
            u = e - nE1
            v = E2[u]
        if v < 0:
            continue
        pu = PHI[u]
        pv = PHI[v]
        if np.isnan(pu) or np.isnan(pv):
            stats[5] += 1
            continue
        dl = int(np.rint(pv - pu))
        if abs(pv - pu - dl) > phi_tol:
            stats[5] += 1
            continue
        stats[0] += 1
        ru, lu = _find(u, parent, pot)
        rv, lv = _find(v, parent, pot)
        if ru == rv:
            if lv == lu + dl:
                stats[2] += 1
            else:
                stats[3] += 1
            continue
        shift = lu + dl - lv            # label rel rv  ->  label rel ru
        ok = _order_ok(S[RAY[v]], S[RAY[v] + 1], v, R, parent, pot, ru, rv, shift, win)
        if ok:
            ok = _order_ok(S[RAY[u]], S[RAY[u] + 1], u, R, parent, pot, ru, rv, shift, win)
        if not ok:
            stats[4] += 1
            continue
        ma = dsum[ru] / dcnt[ru]
        mb = dsum[rv] / dcnt[rv] - shift          # B's mean dev with labels rel ru
        if abs(ma - mb) > comp_tol:
            stats[6] += 1
            continue
        if size[ru] >= size[rv]:
            parent[rv] = ru
            pot[rv] = shift
            size[ru] += size[rv]
            dsum[ru] += dsum[rv] - shift * dcnt[rv]
            dcnt[ru] += dcnt[rv]
        else:
            parent[ru] = rv
            pot[ru] = -shift
            size[rv] += size[ru]
            dsum[rv] += dsum[ru] + shift * dcnt[ru]
            dcnt[rv] += dcnt[ru]
        stats[1] += 1


@njit(cache=True)
def resolve_all(parent, pot, root, lab):
    for x in range(len(parent)):
        r, l = _find(x, parent, pot)
        root[x] = r
        lab[x] = l


@njit(cache=True, parallel=True)
def lis_clean(L, S, nray, keep):
    """Per ray keep the longest strictly increasing subsequence of labels (nodes are in radial order);
    L == INT_MIN means already dropped."""
    for q in prange(nray):
        a0 = S[q]
        a1 = S[q + 1]
        n = a1 - a0
        if n == 0:
            continue
        tails = np.empty(n, np.int64)       # node index of the smallest tail of an increasing run of each length
        prev = np.full(n, -1, np.int64)
        ln = 0
        for i in range(n):
            x = a0 + i
            if L[x] == -2147483648:
                continue
            lo = 0
            hi = ln
            while lo < hi:
                mid = (lo + hi) // 2
                if L[tails[mid]] < L[x]:
                    lo = mid + 1
                else:
                    hi = mid
            if lo > 0:
                prev[i] = tails[lo - 1] - a0
            tails[lo] = x
            if lo == ln:
                ln += 1
        if ln == 0:
            continue
        i = tails[ln - 1] - a0
        while i >= 0:
            keep[a0 + i] = True
            i = prev[i]


@njit(cache=True, parallel=True)
def assign(L, R, S, HS, HR, HW, PIT, nray, reach, D_r):
    """Per ray: for each fitted crossing (radius, winding) the kept sheet node with label == winding, if within reach x
    pitch. D_r[h] = sheet radius - fitted radius (NaN if none). Kept labels increase with radius on every ray."""
    for q in prange(nray):
        a0 = S[q]
        a1 = S[q + 1]
        kl = np.empty(a1 - a0, np.int64)
        kx = np.empty(a1 - a0, np.int64)
        n = 0
        for x in range(a0, a1):
            if L[x] != -2147483648:
                kl[n] = L[x]
                kx[n] = x
                n += 1
        p = PIT[q]
        for h in range(HS[q], HS[q + 1]):
            D_r[h] = np.nan
            if n == 0:
                continue
            i = np.searchsorted(kl[:n], HW[h])
            if i < n and kl[i] == HW[h]:
                d = R[kx[i]] - HR[h]
                if abs(d) <= reach * p:
                    D_r[h] = d


# --------------------------------------------------------------------------------------------- main


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
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
    ap.add_argument('--link-sigma', type=float, default=0.15)
    ap.add_argument('--link-dmax', type=float, default=0.35)
    ap.add_argument('--phi-tol', type=float, default=0.35)
    ap.add_argument('--comp-tol', type=float, default=0.4)
    ap.add_argument('--win', type=int, default=3)
    ap.add_argument('--min-comp', type=int, default=500)
    ap.add_argument('--reach', type=float, default=1.5)
    ap.add_argument('--bbox', type=int, nargs=4, default=None)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--dump', default=None, help='npz of nodes/labels for diagnostics')
    a = ap.parse_args()
    import numba
    numba.set_num_threads(a.threads)
    t0 = time.time()
    S_ = zarr.open(a.surface, mode='r')
    umb = Umbilicus(a.umbilicus, a.umbilicus_scale)
    wind = load_windings(a.windings)
    wids = sorted(wind)
    idx = {w: i for i, w in enumerate(wids)}
    zs = np.arange(a.z0, a.z1, a.dz)
    K, M, NW = len(zs), a.rays, len(wids)
    T = int(a.rmax / a.step)
    tt = np.arange(T, dtype=np.float32) * a.step
    th = np.arange(M) * 2 * np.pi / M
    cos, sin = np.cos(th).astype(np.float32), np.sin(th).astype(np.float32)
    y0, y1, x0, x1 = a.bbox if a.bbox else (0, S_.shape[1], 0, S_.shape[2])
    czc = S_.chunks[0]
    print(f'{NW} windings, {K} slices x {M} rays', flush=True)
    # ---- 1. nodes per slice
    Rl, Pl, Nl, PITl, HRl, HWl, HNl = [], [], [], [], [], [], []
    slab_c, slab = None, None
    kmax = 2 * NW + 16
    for k, z in enumerate(zs):
        c = int(z) // czc
        if c != slab_c:
            slab = None
            slab = np.asarray(S_[c * czc:min((c + 1) * czc, S_.shape[0]), y0:y1, x0:x1])
            slab_c = c
        sl = slab[int(z) - c * czc].astype(np.float32)
        cx, cy = umb.at(float(z))
        mx = (cx - x0 + tt[None, :] * cos[:, None]).astype(np.float32)
        my = (cy - y0 + tt[None, :] * sin[:, None]).astype(np.float32)
        prof = cv2.remap(sl, mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT, borderValue=0)
        cand = np.zeros((M, 1024), np.float64)
        cand_n = np.zeros(M, np.int64)
        sheet_candidates(prof, a.step, a.thr, a.minrun, cand, cand_n)
        hit_r = np.zeros((M, kmax), np.float64)
        hit_w = np.zeros((M, kmax), np.int64)
        hit_n = np.zeros(M, np.int64)
        for w in wids:
            ray_hits(slice_poly(wind[w][1], float(z)), w, cx, cy, M, hit_r, hit_w, hit_n, 120.0)
        cm = int(cand_n.max()) if M else 0
        R = np.zeros((M, max(cm, 1)), np.float32)
        PHI = np.full((M, max(cm, 1)), np.nan, np.float32)
        PIT = np.zeros(M, np.float32)
        HR = np.zeros((M, kmax), np.float32)
        HW = np.zeros((M, kmax), np.int32)
        HN = np.zeros(M, np.int64)
        ray_phi(cand, cand_n, hit_r, hit_w, hit_n, R, PHI, PIT, HR, HW, HN)
        msk = np.arange(R.shape[1])[None, :] < cand_n[:, None]
        Rl.append(R[msk]); Pl.append(PHI[msk]); Nl.append(cand_n.copy()); PITl.append(PIT)
        hm = np.arange(kmax)[None, :] < HN[:, None]
        HRl.append(HR[hm]); HWl.append(HW[hm]); HNl.append(HN.copy())
        if k % 50 == 0:
            print(f'nodes slice {k + 1}/{K} z {z}: {time.time() - t0:.0f} s', flush=True)
    slab = None
    R = np.concatenate(Rl); PHI = np.concatenate(Pl); del Rl, Pl
    cnt = np.concatenate(Nl); S = np.zeros(K * M + 1, np.int64); S[1:] = np.cumsum(cnt)
    PIT = np.concatenate(PITl)
    HR = np.concatenate(HRl); HW = np.concatenate(HWl); del HRl, HWl
    hcnt = np.concatenate(HNl); HS = np.zeros(K * M + 1, np.int64); HS[1:] = np.cumsum(hcnt)
    N = len(R)
    RAY = np.repeat(np.arange(K * M, dtype=np.int64), cnt)
    print(f'{N} sheet nodes, {len(HR)} fitted crossings, phi finite {np.isfinite(PHI).mean():.3f}: '
          f'{time.time() - t0:.0f} s', flush=True)
    # ---- 2. links
    E1 = np.full(N, -1, np.int64); C1 = np.full(N, np.inf, np.float32)
    E2 = np.full(N, -1, np.int64); C2 = np.full(N, np.inf, np.float32)
    build_links(R, S, PIT, K, M, a.link_dmax, a.link_sigma, E1, C1, E2, C2)
    Cc = np.concatenate([C1, C2]); del C1, C2
    nlinks = int(np.isfinite(Cc).sum())
    order = np.argsort(Cc, kind='stable')[:nlinks]
    del Cc
    print(f'{nlinks} links (ray {int((E1 >= 0).sum())}, slice {int((E2 >= 0).sum())}): {time.time() - t0:.0f} s',
          flush=True)
    # ---- 3. union-find
    parent = np.arange(N, dtype=np.int64)
    pot = np.zeros(N, np.int64)
    size = np.ones(N, np.int64)
    st = np.zeros(7, np.int64)
    fin0 = np.isfinite(PHI)
    dsum = np.where(fin0, PHI, 0.0).astype(np.float64)
    dcnt = fin0.astype(np.float64) + 1e-9
    union_all(order, N, E1, E2, PHI, R, RAY, S, parent, pot, size, dsum, dcnt, a.phi_tol, a.comp_tol, a.win, st)
    del order, E1, E2, dsum, dcnt
    ust = dict(zip(['used', 'merged', 'cycle_ok', 'cycle_bad', 'refused_order', 'skipped_phi', 'refused_dev'],
                   map(int, st)))
    print('union', ust, f'{time.time() - t0:.0f} s', flush=True)
    root = np.empty(N, np.int64); lab = np.empty(N, np.int64)
    resolve_all(parent, pot, root, lab)
    del parent, pot
    # ---- 4. one offset per component
    val = PHI - lab
    fin = np.isfinite(val)
    o = np.lexsort((val[fin], root[fin]))
    rr = root[fin][o]; vv = val[fin][o]
    starts = np.r_[0, np.nonzero(np.diff(rr))[0] + 1]
    ends = np.r_[starts[1:], len(rr)]
    med = vv[(starts + ends) // 2]
    base = np.zeros(N, np.int64); has = np.zeros(N, bool)
    base[rr[starts]] = np.rint(med).astype(np.int64); has[rr[starts]] = True
    csize = np.bincount(root, minlength=N)
    good_c = has & (csize >= a.min_comp)
    L = np.where(good_c[root], lab + base[root], -2147483648).astype(np.int64)
    dev = np.abs(val - base[root])
    comp_ids = np.nonzero(good_c)[0]
    big = np.sort(csize[comp_ids])[::-1]
    cst = dict(components_kept=int(len(comp_ids)), nodes_in_kept=int(good_c[root].sum()),
               share_in_kept=round(float(good_c[root].mean()), 4), largest=[int(x) for x in big[:10]],
               fit_disagrees_gt_half=round(float((dev[good_c[root] & fin] > 0.5).mean()), 4))
    print('components', cst, f'{time.time() - t0:.0f} s', flush=True)
    if a.dump:
        np.savez(a.dump, R=R, PHI=PHI, RAY=RAY, root=root, lab=lab, base=base, M=M, K=K, PIT=PIT)
    del val, dev, root, lab
    # ---- 5. radial order per ray
    keep = np.zeros(N, np.bool_)
    lis_clean(L, S, K * M, keep)
    nlab = int((L != -2147483648).sum())
    L = np.where(keep, L, -2147483648)
    print(f'labelled {nlab}, kept after radial order {int(keep.sum())}', flush=True)
    # ---- 6. assign windings
    D_r = np.full(len(HR), np.nan, np.float64)
    assign(L, R, S, HS, HR, HW, PIT, K * M, a.reach, D_r)
    print(f'crossings matched {int(np.isfinite(D_r).sum())} of {len(HR)}: {time.time() - t0:.0f} s', flush=True)
    F = np.full((NW, K, M), np.nan, np.float32)
    hray = np.repeat(np.arange(K * M, dtype=np.int64), hcnt)
    lut = np.full(int(max(wids)) + 2, -1, np.int64)
    for w in wids:
        lut[w] = idx[w]
    wi = lut[np.clip(HW, 0, len(lut) - 1)]
    okh = (wi >= 0) & np.isfinite(D_r)
    F[wi[okh], hray[okh] // M, hray[okh] % M] = D_r[okh]
    # ---- 7. move nodes (matched only; everything else a hole)
    os.makedirs(a.out, exist_ok=True)
    rep = dict(args=vars(a), union=ust, components=cst, nodes=int(N), windings={})
    cxy_z = np.array([umb.at(float(z)) for z in zs])
    for w in wids:
        name, P0, meta = wind[w]
        P = densify(P0, a.densify)
        wj = idx[w]
        H, W, _ = P.shape
        flag = np.full((H, W), 3, np.uint8)
        dout = np.zeros((H, W), np.float32)
        ok = np.isfinite(P).all(-1)
        Z = P[..., 2]
        inz = ok & (Z >= zs[0] - a.dz / 2) & (Z <= zs[-1] + a.dz / 2)
        rr_, cc_ = np.nonzero(inz)
        z = Z[rr_, cc_]
        k = np.clip(np.rint((z - zs[0]) / a.dz).astype(int), 0, K - 1)
        cxy = cxy_z[k]
        dx, dy = P[rr_, cc_, 0] - cxy[:, 0], P[rr_, cc_, 1] - cxy[:, 1]
        ang = np.mod(np.arctan2(dy, dx), 2 * np.pi)
        m = np.rint(ang / (2 * np.pi) * M).astype(int) % M
        d = F[wj, k, m]
        good = np.isfinite(d)
        flag[rr_[good], cc_[good]] = 0
        r = np.hypot(dx, dy)
        scale = np.where(r > 1e-6, (r + np.nan_to_num(d)) / np.maximum(r, 1e-6), 1.0)
        Pn = np.full_like(P, -1.0)
        Pn[rr_[good], cc_[good], 0] = (cxy[:, 0] + dx * scale)[good]
        Pn[rr_[good], cc_[good], 1] = (cxy[:, 1] + dy * scale)[good]
        Pn[rr_[good], cc_[good], 2] = P[rr_[good], cc_[good], 2]
        dout[rr_[good], cc_[good]] = d[good]
        od = os.path.join(a.out, name)
        os.makedirs(od, exist_ok=True)
        for i, c in enumerate('xyz'):
            tifffile.imwrite(os.path.join(od, f'{c}.tif'), Pn[..., i].astype(np.float32))
        tifffile.imwrite(os.path.join(od, 'snap_flag.tif'), flag)
        tifffile.imwrite(os.path.join(od, 'snap_d.tif'), dout)
        meta = dict(meta)
        sc = meta.get('scale', [0.05, 0.05])
        meta['scale'] = [sc[0] * a.densify, sc[1] * a.densify]
        meta['source'] = f"label_snap of {meta.get('source', 'fit')}"
        json.dump(meta, open(os.path.join(od, 'meta.json'), 'w'), indent=1)
        rep['windings'][name] = dict(nodes=int(ok.sum()), in_z=int(inz.sum()), moved=int(good.sum()),
                                     median_abs_move=float(np.median(np.abs(d[good]))) if good.any() else None)
    rep['total'] = {k2: sum(v[k2] for v in rep['windings'].values()) for k2 in ('nodes', 'in_z', 'moved')}
    rep['seconds'] = round(time.time() - t0, 1)
    json.dump(rep, open(os.path.join(a.out, 'report.json'), 'w'), indent=1)
    print('DONE', json.dumps(rep['total']), f'{time.time() - t0:.0f} s', flush=True)


if __name__ == '__main__':
    main()
