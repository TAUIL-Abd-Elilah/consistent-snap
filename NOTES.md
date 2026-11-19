# Consistent snap — pre-registration and results log (times UTC, written as the work happened)

## Why this (18 Nov 2026, ~22:20 UTC)
- Paul, #unrolling-vc3d, 3 Nov 16:26 UTC: snapping is "trickier than it seems, and unsolved"; the hard case is
  snapping a whole scroll consistently so that wraps do not coincide after snapping (paraphrased).
- sean, 23 Oct: First Letters submissions mostly fail on bad or too-small segmentations; a spiral fit gives huge areas of
  decent surface. Prize page: unwrapping must show fibres and no sheet jumps in cross-section; real data.
- villa: no consistent snapper. Chase's per-label snap_labels (#2368) closed by the inactivity bot on 18 Nov; lasagna
  offset_snap broken (#2391, #2418 open).
- Our own 4 Nov plan for this was never executed until now.

## Method (consistent_snap.py)
Per z slice (every dz), 2048 rays from the umbilicus; fitted windings' crossings vs. predicted-sheet runs along each ray;
exact order-preserving one-to-one alignment (DP); two passes (raw, then around a NaN-median prior); nodes moved radially.
`--mode naive` = same pipeline, nearest sheet per winding (the per-surface baseline).

## Known-answer test, fixed BEFORE the first score (18 Nov ~22:50 UTC, clock from date)
- Data: PHerc0139 volume 20260102150214, fit frame L2 (9.596 um), band z 3200-4200.
- Fit: our arm-A seed-1 fit `novA1` (villa fit_spiral, no patches; Nov fit A/B), 115 windings.
- Sheet evidence: public recto prediction `20260102150214-surface-20260203205514-surface-recto-2um-ps192-L0-th0.5`, level 2.
- Truth: the 36 human wrap segments w023-w059 (minus w046), never an input; scored with the unchanged Nov A/B scorer
  `score_band.py` (within 2 / 4 voxels of a winding surface, adjacent-pair exactness, winding hops, per-wrap modal share,
  duplicate-sheet share).
- Parameters fixed (defaults in the script at this time): thr 127, dz 4, rays 2048, step 0.5, densify 4, sigma 0.35,
  dmax 0.6, dmax2 0.45, cs 1.0, median window 3 x 9, gaussian 1.5 x 6. Set from slice views of z 3560-3680 (prediction
  and fit only, no wraps shown).
- Change before any score (18 Nov 22:52 UTC): densify 2 instead of 4 for the scored runs (4x gives ~60 M triangles
  per arm, beyond the free RAM of the scorer). No score had been computed.
- Change before any score (18 Nov 22:57 UTC): the scorer's duplicate share depends on node density (12 nearest
  neighbours within 25 vx); at densify 2 most of those neighbours are the same winding, which would flatter C and N.
  So C and N are also scored at the fit's own node spacing: the snapped grid taken every 2nd node ([::2, ::2]), which
  are exactly the fit's own nodes after the move (131 x 307 -> 66 x 154, corner-aligned). The duplicate bar uses this
  native-spacing C; within_2 and adjacent-exact bars use the densified C as fixed above. Both are reported.
- Arms: A = raw fit, N = naive snap, C = consistent snap.
- Bars (decided now): C within_2 >= A within_2 + 15 points; C duplicate share <= A duplicate share; C adjacent-exact
  >= A adjacent-exact - 2 points. N is reported, not barred. If C fails, the result is reported as failed.

## v2: track_snap.py (sweep), fixed 18 Nov 23:12 UTC, before any v1 or v2 score was seen
- Why: a CT view of v1 (z 3700, crop y 2480-2880 x 4150-4550, fit / v1 / human wraps) and the truth-free step count
  showed that v1 turns the fit's drift across sheets into stair-steps: 8.3 % of in-slice neighbour pairs move > 5 vx
  apart (4.1 % > 8 vx). Sheet jumps in cross-section are what the prize page says an unwrap must not have.
- Method: one joint sweep (middle slice, then up, then down; rays from m0 = 0). Per ray the same exact order-preserving
  one-to-one alignment, but each winding's target is where its sheet continues (z-carry from the previous visited slice,
  else ray-carry along the slice, else its fitted radius), with a tight tolerance; displacements may drift up to 2.5
  pitches from the fit. First try (independent Viterbi per winding, then drop collisions) was abandoned on the smoke
  slab before any score: 3.5 M of 7.0 M cells collided.
- Parameters fixed (script defaults now): sig-a 0.35, dmax-a 0.6, sig-t 0.15, dmax-t 0.35, cs 1.0, carry 16, reach 2.5,
  m0 0, densify 2, med 3 x 9, sig 1.5 x 6, thr 127, dz 4, rays 2048, step 0.5. Set a priori; checked once on the smoke
  slab z 3560-3680 (CT + fit + v1 + v2 at z 3620, no wraps shown).
- Arm S = v2 on novA1, full band z 3200-4200. Bars (decided now, same scorer, duplicate share at native spacing):
  S within_2 >= A + 15 points; S duplicate share <= A; S adjacent-exact >= A - 2 points; and, for sheet switches,
  S hop rate (winding_hops_total / within_4 nodes) <= A hop rate and S mean per-wrap modal share >= A.
- Replication required before any claim: the same parameters on novA2 (seed 2), same bars against its own raw fit.

## Results (18 Nov, 23:40 UTC; human-wrap scores still running)
- v1 C FAILS its duplicate bar. Duplicate share at native spacing (score_ab.duplicate_sheet_share, unchanged):
  A 3.6 %, C 18.2 %, N 26.2 %, C with steps cut (cut_steps.py, 6 vx) 10.6 %, sweep 31.3 %, grow 26.0 %.
  Source: the gap-filled nodes (flag 1, moved by the median field, not one-to-one). C matched nodes only: 3.2 %
  (A on the same nodes 0.4 %); C matched + filled: 24.3 %. Matched nodes are 43 % of C's nodes.
- Truth-free CT check (normal_profile.py, box z 3500-3600 y 3000-3600 x 2200-2800, not looked at before):
  contrast I(0) - I(+-4..6): A -1.6 (flat, off sheet), C +17.6, N +19.5, sweep +18.5, grow +15.7, human wraps +21.6.
  Snapped surfaces sit on the recto face (peak at +1..+2 vx along the normal, air-side trough at -3).
- NEW ARM fixed 18 Nov 23:45 UTC, before any human-wrap score of any snap arm was seen: Cm = C keeping only the
  nodes whose own ray found a one-to-one match in pass 2 (flag 0); everything else becomes a hole (-1). Rationale from
  the truth-free duplicate diagnosis above: a winding left unmatched by the alignment usually lost its sheet to a
  neighbour (the fit has more windings than sheets there), so filling it re-creates the stacking the method exists to
  prevent. Bars for Cm (post hoc, so replication is required: novA2 and band 2 b2A1 with the same rule):
  Cm duplicate share (native) <= A; Cm within_2 >= A within_2 (no recall loss despite the holes); Cm adjacent-exact
  >= A - 2 points; and precision (share of Cm nodes within 2 vx of a human wrap, among nodes within 12 vx of one)
  >= A precision + 15 points.
- Cm on novA1, human wraps (score_band.py unchanged, B-only run, finished 19 Nov ~02:20 UTC; paired bootstrap against
  arm A rows of score/ab_seed1.json with score_ab.paired_bootstrap): within_2 31.3 -> 52.2 % (+20.9 points, 95 % CI
  +18.4 to +23.5); within_4 60.4 -> 72.8 %; adjacent-exact 63.9 -> 66.8 %; same-winding 17.9 -> 17.9 %; hop rate
  (hops / within_4 nodes) 0.0272 -> 0.0206; mean per-wrap modal share 0.161 -> 0.159 (unchanged: the stair-steps
  still switch winding labels along a wrap); duplicate share at native spacing 3.6 -> 3.2 % (scorer's own value on the
  densified mesh 4.7 %, density-inflated as noted above). Bars: within_2 PASS, duplicate PASS, adjacent PASS,
  precision PASS (tools/precision.py, 400k sampled nodes, zone 12 vx: A 28.9 % within 2 vx, Cm 67.4 %, C 58.8 %,
  N 62.7 %; within 4 vx: A 54.5, Cm 87.0, C 79.4, N 82.5). Replication duplicate shares (native): novA2 3.5 -> 3.2 % PASS; band 2 b2A1 2.29 -> 2.44 % FAIL
  by 0.15 points.
- Jumps (second difference of the move > 10 vx between densified neighbours): C 6.4 % in-slice / 6.0 % z,
  N 8.9 / 8.3, sweep 12.5 / 7.4, grow 9.0 / 8.7. The two tracking variants are worse than v1 on every truth-free
  measure: dropped (code kept in track_snap.py).
- Replication precision (tools/precision.py, same settings; finished 19 Nov ~04:35 UTC): novA2 A 28.8 % -> Cm 67.4 %
  within 2 vx (within 4: 54.4 -> 87.0); band 2 b2A1 A 27.7 % -> Cm 67.1 % (within 4: 52.8 -> 87.1). Precision bar
  (>= A + 15 points) PASS on both. Human-wrap recall scores for Cm_A2 / Cm_b2 still running (started 02:22 UTC).
- Replication human-wrap scores (score_band.py / score_band_b2.py unchanged, B-only; compare_arms.py vs ab_seed2 /
  ab_b2 arm A; finished 19 Nov ~06:40 UTC):
  novA2: within_2 31.5 -> 52.9 % (+21.4, CI +18.9 to +24.1); within_4 61.1 -> 73.5 %; adjacent-exact 64.0 -> 66.9 %;
  hop rate 0.0269 -> 0.0205; modal 0.167 -> 0.165. Bars: within_2 PASS, adjacent PASS, duplicate PASS, precision PASS.
  b2A1: within_2 30.6 -> 53.4 % (+22.8, CI +19.7 to +26.1); within_4 60.3 -> 74.8 %; adjacent-exact 67.3 -> 67.6 %;
  hop rate 0.0148 -> 0.0206 (WORSE; not a Cm bar, reported); modal 0.117 -> 0.127. Bars: within_2 PASS, adjacent
  PASS, precision PASS, duplicate FAIL (2.29 -> 2.44 %, native).
- Verdict under the pre-set rule: Cm replicates on recall, precision and adjacency on all three fits; the duplicate
  bar holds on two of three (band 2 misses by 0.15 points). Stair-steps remain (modal share flat; band-2 hop rate up).
- Pre-registered arms C and N on novA1 (score_band.py, finished 19 Nov 06:33-08:08 UTC). C (densified): within_2
  31.3 -> 64.1 % (+32.8, CI +30.6 to +35.2), within_4 85.4 %, adjacent-exact 67.1 %, hops 520,508, duplicate 17.1 %
  (native 18.2 %) -> C FAILS its pre-set duplicate bar, as already known from the truth-free check. C native spacing:
  within_2 61.8 %, adjacent 67.4 %. N (naive, densified): within_2 63.2 %, adjacent 65.4 %, hops 1,102,441, duplicate
  25.7 %. So consistent vs naive: same recall, about half the hops and a third fewer duplicates; Cm trades 12 points of
  recall (64 -> 52 %) for no duplicates.
  N native spacing: within_2 61.8 %, adjacent 66.1 %, hops 735,211, duplicate 26.2 % (all scoring runs now finished).

- Published 19 Nov ~10:57 UTC: this repo, commit 1c7e5d9, release v1 (fits + C/Cm meshes).

## v3: label_snap.py (sheet-graph winding labels), bars fixed 19 Nov 2026 14:19 UTC (date -u), before any code was run
- Why: Cm passes, but each human wrap still switches fitted winding as often as with the raw fit (modal share flat;
  band-2 hop rate up). Cause: the per-ray alignment picks the sheet nearest the fitted winding, so where the fit drifts
  across sheets the snapped winding steps to the next sheet.
- Method: sheet nodes = predicted-sheet runs on every ray of every slice (as in v1). Link nodes of neighbouring rays
  and neighbouring slices by the same exact order-preserving one-to-one alignment (tight tolerance). Union-find with
  integer label potentials, best links first; the label jump across a link is round(phi(v) - phi(u)), phi = the fit's
  fractional winding id interpolated along the ray (0 except at the fit's winding seam). A link is refused if it would
  break the radial order of labels among nearby sheets on either ray. Each sheet component gets ONE integer offset to
  the fit (median of phi - potential), so a sheet keeps one winding label however far the fit drifts. Per ray,
  labels that are not strictly increasing outward are dropped. Fitted winding w moves to the sheet labelled w on its
  ray if within a reach of 1.5 local pitches; everything else is a hole (as in Cm).
- Parameters are set on the smoke slab z 3560-3680 from CT/prediction/fit views and truth-free counts only (no
  wraps shown), then frozen in this file before the first human-wrap score.
- Arm U = label_snap on novA1 (z 3200-4200), densify 2; duplicate share at native spacing, same scorer and tools.
  Bars (each against the same fit's raw arm A and its Cm):
  U mean per-wrap modal share >= Cm + 0.03; U hop rate <= 0.75 x Cm hop rate;
  U within_2 >= A + 15 points; U duplicate share (native) <= A; U adjacent-exact >= A - 2 points.
  Replication with the same frozen parameters on novA2 and b2A1 is required before any claim; a fit that fails any
  bar is reported as a failure.
- Development (smoke slab z 3560-3680, truth-free only; no wraps shown or scored), 19 Nov 14:20-14:27 UTC:
  first run without a component check: components crept across sheets (a component spanning radius 319-1116 with
  labels -2..3, phi - label spread +-6), only 14 % of fit nodes moved. Added --comp-tol: a merge is refused if the two
  components' mean (phi - label) differ by more than comp-tol after the label shift. With 0.4: dev spread within the
  largest components +-0.5, 32 % of in-band fit nodes moved (Cm on the same slab: 44 %), jump share (2nd difference
  > 10 vx) in-slice 3.2 % / z 2.3 % (Cm 4.7 / 4.2), duplicate share native 0.66 % (Cm 2.3 %). CT view
  (local view, not in the repo): U lines follow sheets with fewer zigzags; steps remain in the dense lower region.
  Sweep (phi-tol, comp-tol) = (0.2,0.4) (0.35,0.6) (0.2,0.6) (0.35,0.3): more coverage costs more jumps; rule set
  before the sweep: most moved nodes with jump share <= 3.5 % both axes and duplicate <= 1 % -> (0.35, 0.4).
- FROZEN 19 Nov 14:27 UTC (date -u), before any human-wrap score of U: label_snap.py defaults (dz 4, rays 2048, step 0.5,
  thr 127, densify 2, link-sigma 0.15, link-dmax 0.35, phi-tol 0.35, comp-tol 0.4, win 3, min-comp 500, reach 1.5),
  same bbox and umbilicus per band as v1. Bars as fixed at 14:19 UTC above.
- U runs (19 Nov ~14:30 UTC, ~100 s each): moved 2.76 M / 2.80 M / 2.82 M of 7.84 / 8.11 / 8.14 M in-band nodes
  (35 %; Cm 43 %). Truth-free: duplicate share at native spacing U 1.19 % (novA1), 1.16 % (novA2), 1.08 % (b2A1),
  all below A (3.6 / 3.5 / 2.29 %) -> duplicate bar PASS on all three (band 2 included, where Cm failed). Area
  186 / 189 / 183 cm2 (scorer's measure). Jump share novA1 full band: U 5.2 % in-slice / 3.8 % z, Cm 5.6 / 4.8.
  Human-wrap scores running (started ~14:32 UTC).

### v3 U human-wrap scores (scorers finished ~16:18 UTC 19 Nov 2026, date -u)

| | seed 1 | seed 2 | band 2 |
|---|---|---|---|
| within 2, A -> U (CI) | 31.3 -> 45.2 (+13.8, +11.2..+16.6) | 31.5 -> 46.5 (+14.97, +12.4..+17.8) | 30.6 -> 45.8 (+15.1, +12.4..+18.1) |
| within 4 | 60.4 -> 62.7 | 61.1 -> 64.5 | 60.3 -> 64.4 |
| adjacent exact | 63.9 -> 58.0 | 64.0 -> 59.3 | 67.3 -> 57.3 |
| same-winding pairs | 17.9 -> 20.3 | 17.8 -> 19.9 | 17.8 -> 21.2 |
| hop rate (Cm) | 0.0109 (0.0206) | 0.0109 (0.0205) | 0.0123 (0.0206) |
| modal share (Cm) | 0.178 (0.159) | 0.178 (0.165) | 0.147 (0.127) |
| duplicate native vs A | 1.19 vs 3.6 | 1.16 vs 3.5 | 1.08 vs 2.29 |

Bars: modal >= Cm + 0.03 FAIL x3 (+0.019, +0.013, +0.020); hop rate <= 0.75 Cm PASS x3 (0.53, 0.53, 0.60 of Cm);
within_2 >= A + 15 FAIL seed 1 (+13.8), FAIL seed 2 (+14.97), pass band 2 (+15.1); duplicate PASS x3;
adjacent >= A - 2 FAIL x3 (-5.9, -4.7, -10.0).
**Verdict: U FAILS its pre-registered bars. Kept in experiments/ as a reported failure; Cm (v1) stays the released arm.**
Reading: winding hops halve and stacking drops to ~1 %, so the labels do hold sheets together, but whole components
carry a wrong integer offset (base = round(median(phi - label))): neighbouring wraps then land on the same or
non-adjacent windings (same-winding pairs up, adjacent exact down 5-10 points). The stair-step became an offset error
between components, not a fix.
