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

## v4: surface renders, raw fit (A) vs matched-only snap (Cm), rules fixed 19 Nov 2026 20:48 UTC (date -u), before any render was made
- Purpose: show what the snap changes for reading: CT rendered on the surface. No new method, no parameter tuning.
- Boxes (the figure crops already published, fixed before): seed 1 and seed 2 z 3400-4000, y 2950-3350, x 4250-4650;
  band 2 z 4700-5100, y 2860-3260, x 4220-4620.
- Pairing: Cm meshes are the raw fit densified x2, moved radially node by node, so grid index = same fit position.
  The raw fit is densified x2 the same way (consistent_snap.densify). Each winding with >= 2000 Cm-valid nodes in the
  box is rendered for both arms on the SAME pixels: pixels where the Cm node is valid (holes masked for both arms).
- Render: each densified grid upsampled bilinearly to ~1 CT voxel per pixel; normal from the grid; CT level 1
  (4.8 um) sampled trilinearly along the normal at offsets -3..+3 level-1 voxels; images: offset-0 slice and the
  max over the 7 offsets. Display: same fixed window for both arms (box 1st-99th percentile of CT).
- Truth-free numbers on the paired pixels, per winding and pooled (not bars; reported as is): mean offset-0
  intensity; share of pixels whose offset-0 intensity is above the box median; texture = std of the
  Laplacian of the offset-0 image. Expect Cm > A on all (consistent with the normal-profile check); any winding where
  A >= Cm is listed.
- Human wraps rendered in the same boxes as a reference only (separate parametrization, not paired).
- Which windings go in the figure: the 3 with the most paired pixels per box, no picking by look.
- Amendment 20:51 UTC (date -u), before any fit render existed (first run rendered only the human wraps,
  which are reference only): the >= 2000-node rule left 0 windings (max 1974 in box 1, 1215 in band 2; the box is
  small). New rule: windings with >= 300 Cm-valid nodes in the box. Everything else unchanged.
- Results 21:00 UTC (date -u), all 3 fits (b2 rerun with its own wrap folder, the band-2 scorer's; fit
  renders unchanged, pooled numbers identical). Paired pixels, offset-0 CT level 1 (A -> Cm):
  seed 1: 30 windings, 15.6 M px; mean 78.3 -> 96.3; above box median 59.7 -> 82.5 %; texture 11.6 -> 14.5.
  seed 2: 30 windings, 15.4 M px; mean 78.5 -> 96.7; above median 59.9 -> 82.7 %; texture 11.9 -> 14.9.
  band 2: 23 windings, 10.4 M px; mean 70.6 -> 94.9; above median 55.2 -> 85.2 %; texture 10.5 -> 13.8.
  Windings where A >= Cm on any of the three: none (0 of 83).
  Human wraps in the same boxes (reference): mean 86.3 / 78.5, above median 70.4 / 64.7 %, texture 12.6 / 12.5.
  Cm is brighter than the wraps: it sits nearer the sheet's dense core, the wraps nearer the face. Not a claim that
  Cm beats human tracing; ink is not tested here. Visible in the renders: papyrus fibre texture where the raw fit
  shows dark gaps between sheets; dark seams in Cm are the remaining stair-steps.
  Figures figures/render_{seed1,seed2,band2}.jpg: top 3 windings by paired pixels, plus the wrap with most pixels in the box.

## v5: label_snap --stitch (joint offsets between components), bars fixed 19 Nov 2026 21:05 UTC (date -u), before any code was written
- Why: U failed because each sheet component took its own offset round(median(phi - label)); neighbouring pieces
  got offsets off by one, so steps moved between pieces (adjacent exact -5..-10, same-winding pairs up).
- Method: everything in U unchanged up to the components. New step: on every ray, two radially consecutive sheet
  nodes from different kept components a, b, with a radial gap below gap-max x local pitch (no missing sheet
  between), vote that b's offset = a's offset + (1 - lab_b + lab_a), i.e. consecutive sheets on a ray are consecutive
  windings. Per component pair, the modal vote is used if it has >= vote-min votes and >= vote-share of the pair's
  votes. Components are joined greedily, most votes first, by union-find with integer potentials; a join that
  contradicts an earlier one is refused. Each joined group gets ONE offset to the fit: round(median(phi - label))
  over all its nodes. Then LIS and assignment exactly as in U.
- Parameters (gap-max, vote-min, vote-share) set on the smoke slab z 3560-3680 from truth-free counts only (share of
  consecutive kept node pairs whose labels differ by exactly 1; share of nodes with |phi - L| > 0.5; nodes moved;
  jump shares; duplicate share), then frozen here before the first human-wrap score. U's frozen parameters stay.
- Bars (same five as v3, against the same fit's A and Cm), on ALL THREE fits (novA1, novA2, b2A1); any failure is
  reported as a failure:
  modal share >= Cm + 0.03; hop rate <= 0.75 x Cm hop rate; within_2 >= A + 15 points; duplicate share (native)
  <= A; adjacent-exact >= A - 2 points.
- Development 21:07 UTC (smoke slab, truth-free): U reproduces (moved 301,043, as before). Stitch without a guard:
  consecutive pairs with label diff 1: 0.670 -> 0.870, but one group takes 96 % of nodes, cycle contradictions 1013
  vs 577 consistent, and nodes with |phi - L| > 0.5 go 17 % -> 85 %: votes chain across missed/extra sheets and the
  labels drift away from the fit. Added --group-tol: a join is refused if the two groups' mean (phi - L) differ by more
  than group-tol after the shift (the comp-tol idea one level up).
- Sweep group-tol (smoke, truth-free): 0.25 / 0.4 / 0.6 -> consecutive diff-1 0.701 / 0.711 / 0.754; |phi - L| > 0.5
  18.3 / 20.7 / 31.4 %; moved 304,817 / 302,911 / 296,478 (U 301,043; consecutive 0.670; disagree 17.2 %).
  Diagnosis (U, smoke): radially consecutive kept sheets get label diff 0 / 1 / 2 in 12 / 68 / 20 % of pairs; the fit's
  phi step across those pairs has median 0.69 / 1.03 / 1.38: label skips sit where the fit packs more windings than
  sheets (see pitch over-packing). Labels tied to the fit inherit it; stitching against the fit makes labels drift.
- FROZEN 21:09 UTC (date -u), rule set after seeing the three sweep values above (stated plainly): the largest
  consecutive diff-1 share whose fit disagreement stays <= 25 % -> --stitch --gap-max 1.6 --vote-min 20 --vote-share 0.7
  --group-tol 0.4; all U parameters unchanged.
  Prediction written before scoring: the change from U is small (consecutive 0.67 -> 0.71), so S will very likely FAIL
  the adjacent bar (U was -4.7 to -10 points) and the modal bar. Run anyway because it is cheap and pre-registered.
- S full runs 21:19 UTC: moved 2.76 / 2.82 / 2.84 M (U 2.76 / 2.80 / 2.82); one large group per band (23-28 M
  nodes), cycle contradictions 3801 / 3409 / 3532 vs consistent 1720 / 1597 / 1747; consecutive diff-1 share
  0.634 / 0.649 / 0.635. Duplicate share native S 1.14 / 1.13 / 1.05 % vs A 3.57 / 3.54 / 2.29 % -> PASS x3.
  Human-wrap scorers started ~21:20 UTC.
- 22:00 UTC: moved D:/Snap/runs/* (except *_S) to C:/VesuviusTemp/snap_moved/runs and deleted D:/Snap/ctcache: First Letters' render job needs D: >= 61 GB free and was skipping windings.
- 22:12 UTC: moved D:/Snap/p4 to C:/VesuviusTemp/p4 (D: floor of First Letters' job).
- 20 Nov 00:18 UTC: PC restarted ~23:xx; S scorers and the Paris 4 ink render had died unfinished. Snap outputs and the 11.5 GB ink mirror survived. Restarted both detached (scorers 1 thread each, render 5 threads).
- 05:24 UTC: Paris 4 renders done (paired px 71.5 M; mean ink A 25.6, Cm 28.6, secondary). Normalisation amended to p95 over paired pixels (PREREG_P4 amendment, before any metric). Frame check: our A render correlates 0.22 with the team's published gen-19 ink render (unflipped), 0.01 flipped. Metric job queued last on the First Letters GPU queue (bash C:/VesuviusTemp/p4/metric_job.sh).

### v5 S human-wrap scores (finished ~06:26 UTC 20 Nov, date -u)
| | seed 1 | seed 2 | band 2 |
|---|---|---|---|
| within 2, A -> S (CI) | 31.3 -> 45.2 (+13.8, +11.2..+16.6) | 31.5 -> 46.6 (+15.1, +12.5..+17.9) | 30.6 -> 46.1 (+15.4, +12.7..+18.3) |
| adjacent exact (U) | 63.9 -> 59.8 (58.0) | 64.0 -> 60.1 (59.3) | 67.3 -> 59.0 (57.3) |
| same-winding pairs | 17.9 -> 19.8 | 17.8 -> 19.7 | 17.8 -> 20.5 |
| hop rate (Cm) | 0.0106 (0.0206) | 0.0108 (0.0205) | 0.0117 (0.0206) |
| modal (Cm) | 0.174 (0.159) | 0.176 (0.165) | 0.153 (0.127) |
| duplicate native vs A | 1.14 vs 3.57 | 1.13 vs 3.54 | 1.05 vs 2.29 |
Bars: modal >= Cm + 0.03 FAIL x3 (+0.015, +0.011, +0.026); hop rate PASS x3; within_2 >= A + 15 FAIL seed 1, pass
seeds 2 / band 2; duplicate PASS x3; adjacent >= A - 2 FAIL x3 (-4.1, -3.9, -8.3).
**Verdict: S FAILS, as predicted before scoring.** Stitching moved adjacent exact up 0.8-1.8 points over U, far short.
Stair-step work stops here (decision written 19 Nov: no further tuning).
