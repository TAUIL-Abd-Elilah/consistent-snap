# Paris 4 ink test: does snapping the team's gen-19 spiral fit show more ink?

Rules fixed 19 Nov 2026 21:44 UTC (date -u), before any Paris 4 data was downloaded or rendered.
(Only directory listings, array headers, FIT_NOTES.md/WINDINGS.txt/meta.json of the fit were read.)

## Data (all public, frame = volume 20260411134726, level 2 = 9.6 um)
- Fit A: PHercParis4/segments/20260926044226-spiral-output-generation-19, mesh/spiral-fit (x/y/z.tif, windings.tif;
  coordinates already in level-2 voxels, grid step 20 voxels). Split into per-winding grids by windings.tif
  (ambiguous -2 and invalid -1 nodes dropped).
- Surface evidence: predictions/surfaces/20260411134726-surface-20260413141734-surface-recto-2um-ps256-L0-th0.45.zarr level 2.
- Umbilicus: the spiral dataset umbilicus.json (level-2 coordinates, the fit's own input; the representations/umbilicus file is the same curve at level 0, x4).
- Ink: predictions/ink-3d/20260411134726-ink3d-20260926084000-student-v12A-25000-c7766be9-L1.zarr level 1 (= level 2
  of the scan), the volume the team used for the published gen-19 ink render.

## Band (by rule)
Middle 1536 slices of the fit's z range (4005-17747, centre 10876), aligned to 256-voxel ink chunks: z 9984-11520.

## Arms
- A: gen-19 as published (densified x2 so it shares Cm's grid).
- Cm: consistent_snap.py matched-only, v1 frozen parameters (dz 4, rays 2048, densify 2, all defaults), then mask flag 0.
No parameter is changed for Paris 4.

## Render (emulates spiral-fitting/render_ink.py: vc_render_tifxyz --group-idx 1 --scale 0.25 --num-slices 5, max composite,
then /p95 -> uint8; our own CPU sampler because vc_render_tifxyz hangs on Windows, villa #2562)
- Each winding grid upsampled x5 (1 pixel = 4 level-2 voxels along the surface, as scale 0.25), normal from the grid,
  ink sampled trilinearly at offsets -2,-1,0,+1,+2 level-2 voxels, max.
- Paired pixels only: pixels valid in Cm; A is masked to the same pixels.
- Strip: per-winding renders concatenated left to right in winding order (rows = z), 4-px gaps; each arm normalised by
  its own p95 as render_ink does.

## Metric: the team's get_ink_metrics.py (scrollprize/ink-coverage-32um ensemble, default threshold 0.5), unchanged.

## Bars (Cm against A, same pixels)
1. ink foreground area: Cm >= 1.15 x A;
2. line_score: Cm >= A - 0.05;
3. visual: a fixed crop (the 8192-px slice with the most A ink by the metric's own slice table) shown A | Cm.
Secondary, reported only: mean max-ink on paired pixels; per-slice area table.
If bar 1 or 2 fails, the result is reported as a failure.

## Amendment 20 Nov 2026 05:22 UTC (date -u), before any ink metric was run
The /p95 normalisation was taken over the whole strip, ~72 % of which is empty here (gaps between winding blocks, Cm
holes, band edges), so the 95th percentile sat inside the surface distribution and the strips came out saturated (one
crop looked at: both arms mostly white). render_ink.py's strips are nearly all surface, so its p95 is in effect a
percentile of surface pixels. New rule: p95 over the paired pixels, per arm. Raw renders unchanged. Seen so far:
render_stats.json (mean ink on paired pixels A 25.6, Cm 28.6) and that one crop.

## Technical fix (2026-11-20 19:18 UTC), no change to any rule or bar
The first metric run (19:15 UTC) failed on both arms before predicting: nnU-Net's trainer-class search imports every trainer
module in order and crashed on a legacy v1 file (variants/loss/from_nnunetv1) that needs the absent 'nnunet' package. The
checkpoint needs nnUNetTrainer_250epochs. Our copy of get_ink_metrics.py now skips modules that fail to import (original kept as
get_ink_metrics.orig.py). Same model, weights, folds, inputs and metric. No metric output was produced, so nothing was seen.
Failed logs: failed_1120/. Re-queued at the GPU queue head at 19:18 UTC, after the First Letters band-10 sweep.

## RESULT (2026-11-21 10:50 UTC), scored against the bars above, nothing changed after the first look
Metric run 21 Nov (both arms, 3 folds, same paired pixels, strip 881 x 284,351 px, 71.5 M paired px).
| | A (gen-19 as published) | Cm (snap) | bar | verdict |
|---|---|---|---|---|
| ink foreground | 2,431,652 px (35.86 cm2) | 2,919,156 px (43.04 cm2) | Cm >= 1.15 x A | 1.200 x, PASS |
| line_score | 0.419 | 0.413 | Cm >= A - 0.05 = 0.369 | PASS |
| visual crop (slice x 262144-270336, most A ink) | 426,121 px | 548,403 px | shown A / Cm | shown, letters in both |
Secondary, reported only: column score 0.092 -> 0.128; mean max-ink on paired pixels 25.6 -> 28.6; per-slice table:
of the 30 slices with ink (>20k px in either arm) Cm has more ink in 19 and less in 11, so the gain is not uniform.
Limits: one band (z 9984-11520) of one scroll; our CPU sampler stands in for vc_render_tifxyz; the metric counts
predicted ink area, not legibility. Same text is visible in both arms in the crop; Cm strokes are somewhat fuller.
