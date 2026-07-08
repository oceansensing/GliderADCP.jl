# M38 validation notes — GliderADCP.jl vs the prior Python processing

> 2026-07-07. Full-mission run of `examples/m38_currents.jl` (sound-speed correction from
> payload CTD, QC, IGRF declination, DAC + bottom-track inverse and shear solutions,
> 127 solved yos — matching the 126 yos in the prior processing; the ADCP was
> duty-cycled Nov 3–27 within the Nov–Mar mission).

## Headline internal-quality metrics (GliderADCP.jl inverse)

| Check | Result |
|---|---|
| Dive vs climb consistency (independent half-yo inversions, same DAC) | r_u = 0.983, r_v = 0.986, med \|Δ\| ≈ 2 cm/s (n = 1514 bins) |
| DAC closure | median 0.005 m/s over 127 yos |
| Glider velocity vs unseen bottom track (DAC-only inverse) | r_v = 0.97, med \|Δ\| ≈ 7 cm/s (n = 807) |
| BT-anchored absolute deep water velocity vs inverse bins | med \|Δ\| = 1.6–1.7 cm/s (n = 136k) |
| Shallow bins (z < 30 m) vs surface GPS drift | med \|Δ\| ≈ 4 cm/s (n = 126 yos) |

## The vertical-structure question (and its resolution)

Our sections differ visibly from the prior processing's figures below ~200 m:
early mission (Nov 4–9) we show **subsurface-intensified** (±0.4–0.5 m/s at 250–500 m)
slope-eddy structure; late mission (deep basin) we show near-barotropic per-yo columns,
while the prior figures show smooth **surface-intensified** profiles throughout.

A pure baroclinic sign flip preserves depth means, dive/climb consistency, and DAC
closure — so those checks cannot arbitrate. Three independent arbiters were run:

1. **End-to-end synthetic with depth-varying flow** through the full beam forward model
   (now a permanent regression test): relative velocities exact at every offset;
   inverse recovers du/dz upright and at the right magnitude. Our chain cannot flip
   structure.
2. **Raw transform-level tilt check** (no solver, no reference): binned mean relative
   velocity E_rel(z) within single yos. Since glider velocity is ~constant on average
   over a yo, the *shape* of E_rel(z) ≈ u_ocean(z) + const:
   - Nov 5 / Nov 7 yos: raw tilt (surface − 300 m) = **−0.16 / −0.18 m/s** —
     subsurface intensification is IN THE RAW DATA (real slope-current/eddy signal);
     the prior processing agrees in sign here (−0.06 / −0.18).
   - Nov 22 yo: raw E_rel is **flat** (−0.32 ± 0.01 from 12–760 m) — barotropic.
     Our profile is correspondingly flat; the prior profile imposes a +0.13 m/s
     surface intensification that is **not present in the raw data**.
3. **BT-anchored absolute velocities** (u_rel + u_glider_over_ground from bottom track,
   no DAC/inverse involved) agree with our deep inverse bins to 1.6 cm/s median.

**Conclusion:** GliderADCP.jl sections are faithful to the raw measurements. The prior
Python profiles agree where the signal is strong (early-mission shelf/slope yos — the
matched-yo correlation there is r ≈ 0.8) and diverge where profiles are weak, where
that pipeline's over-smoothing dominates: it ran wSmoothness = 1 at dz = 1 m —
a ~100× stiffer curvature penalty than at the documented dz = 10 m — plus the
documented v2.0.0 transform defects (dive-cast matrix misalignment, halved X/Z,
possible rotation transpose). Mission-wide per-yo correlation against that reference
is therefore low (median r_u ≈ 0.1) *by expectation*, and is not a defect indicator
for this package.

## Notes / future improvements

- Shear-vs-inverse intercomparison: r_u = 0.58, rms ≈ 0.2 m/s pooled. The shear path
  currently uses simple depth-mean DAC referencing; per-yo reference offsets and
  integration drift dominate the discrepancy. Planned: time-in-bin-weighted referencing
  (gliderad2cp semantics) and per-cast integration.
- QC rejects 52 % of beam samples on the full mission with default thresholds (SNR
  floor + amplitude + correlation + velocity cap + surface mask); revisit per-screen
  contributions when tuning for deep, quiet water.
- Surface-drift comparison is noisy (windage/Stokes on a surfaced glider); it bounds
  gross errors but should not be over-interpreted below ~5 cm/s.

## Why the shear and inverse methods disagree (2026-07-07 diagnosis)

Pooled shear-vs-inverse agreement on M38 is r_u ≈ 0.58, rms ≈ 0.22 m/s, worst near the
surface. Quantitative diagnosis (script results, full mission):

1. **A measured range-dependent velocity bias.** The mission-mean *within-ping* relative
   velocity, rotated into the glider track frame, decays with range: along-track slope
   **−3.3×10⁻⁴ s⁻¹ (dives) and −3.1×10⁻⁴ s⁻¹ (climbs)** — same sign and size for both
   cast directions — while the cross-track component is 10–100× smaller (−2×10⁻⁵ / +3×10⁻⁶).
   Per cell the bias is tiny (±3.5 mm/s across the 30-m window); it is invisible to any
   per-sample QC. This is the classic glider-ADCP "shear bias" (Todd et al. 2017;
   gliderad2cp's `process_bias`, whose velocity-dependent variant exists for exactly this
   speed-correlated, track-aligned signature).
2. **The shear method integrates it; the inverse does not.** Integrating −3.2×10⁻⁴ s⁻¹
   over a 500-m yo predicts a −0.16 m/s top-to-bottom tilt; the observed per-yo tilt of
   (u_shear − u_inverse) is **median −4.6×10⁻⁴ s⁻¹ ≈ −0.23 m/s per 500 m** — same sign,
   same order. DAC referencing then pivots the tilted profile about the yo's
   (time-weighted) mean depth, so the shear solution comes out **+0.13 m/s too fast near
   the surface and −0.10 m/s too slow at depth** (median u_inv − u_sh: −0.149 at 0–50 m,
   +0.105 at 600–1000 m). In the inverse, each sample constrains only its own depth bin
   against the per-ping glider unknown — the bias adds ~mm/s per bin instead of
   integrating.
3. **Independent surface arbiter.** Against surface GPS drift (z < 30 m bins):
   inverse med|Δ| = 0.042 m/s, bias −0.025; **shear med|Δ| = 0.151 m/s, bias +0.111** —
   the shear top is confirmed biased, the inverse is not.
4. **Integrated noise.** Beyond the tilt, integrating per-bin shear noise random-walks:
   rms(u_inv − u_sh) is 0.18–0.25 m/s at all depths, vs ~2 cm/s per-bin consistency for
   the inverse (dive/climb check). The upper column additionally has ~3× fewer samples
   per bin (median 46–67 vs ~180 at mid-depth) — sampled only during dive starts and
   climb ends.
5. **Inverse outliers noted:** 289 bins reach |u| up to 1.6 m/s at 345–605 m, all in
   yos 8–18 (the slope-eddy transit, DAC up to 0.51 m/s). Cast-consistent (dive/climb
   r = 0.98) so not a transform artifact, but flagged for the Phase-7 compass-deviation
   check.

**Conclusion:** the disagreement is a property of the *shear method* inheriting an
integrable instrument bias plus integrated noise — not of the common trunk (shared by
both methods) nor the inverse. The inverse (+ bottom track) sections are the reference
product. **Phase-7 item 1:** shear-bias correction — subtract the mission-calibrated
along-track bias profile b(r) before differencing (real ocean shear averages out of the
track-frame mission mean under varied headings), and/or the gliderad2cp
displacement-regression form; acceptance = shear-vs-inverse tilt collapse and drift-test
parity.

## Shear-bias correction: implementation and what it does (and doesn't) fix (2026-07-07)

Implemented as `shear_bias` / `apply_shear_bias!` / `calibrate_shear_bias!`:

- **Estimator**: mission-mean of the *adjacent-pair differences* of track-frame relative
  velocities — the exact sample population the shear estimator consumes (with partial
  coverage, mean-of-differences ≠ difference-of-means, so the earlier per-offset mean
  profile under-corrected). Pair differences are integrated to a mean-removed bias
  profile B(offset). An optional `velocity_scaled=true` form fits bias ∝ ping
  through-water speed (it did not outperform the plain form on M38).
- **Application** re-demeans the subtracted profile over each ping's finite cells, so
  ping-mean velocities are unchanged **by construction** — the inverse's glider-velocity
  and DAC content is untouched (verified: median profile change 3 mm/s).
- **M38 calibration**: slope −3.9…−4.1×10⁻⁴ s⁻¹ along-track, cross-track ~10⁻⁶,
  heading concentration R = 0.23 (safely diverse). One pass removes it to machine zero,
  and the **residual pairwise bias is < 1×10⁻⁴ s⁻¹ in every cell-depth band** (was
  ≈ −4×10⁻⁴ uniformly; also verified depth-stratified: −3.0…−3.8×10⁻⁴ before).

Effect on M38 solutions: shear-vs-inverse r_u 0.582 → 0.648, rms 0.217 → 0.197 m/s;
per-yo tilt median −4.6 → −3.3×10⁻⁴ s⁻¹; surface-drift bias of the shear top
+0.111 → +0.087 m/s. Dive/climb and DAC-closure metrics of the inverse unchanged.

**What the remaining disagreement is (investigated to closure):** after correction the
shear *samples* are provably unbiased (< 10⁻⁴ s⁻¹ at every depth), the bin statistic is
not the cause (mean and median give identical tilt), and the discrepancy vs the inverse
concentrates at 200–700 m as per-yo *integration drifts* of ±0.1–0.3 m/s (e.g. a yo
whose raw relative velocity is constant to ±0.01 from 215–505 m shows a −0.2 m/s shear
drift while the inverse stays flat). The residual median tilt (−3.3×10⁻⁴ s⁻¹ ≈ 3× the
median's standard error given the ±1.3×10⁻³ IQR) is the noise-skew of those random
walks, not a removable sample bias. This is the intrinsic error-propagation difference
between the methods: **integration accumulates what the inverse localizes.** The shear
product is retained as a corrected, second-opinion diagnostic; the inverse (+ bottom
track) remains the reference. Possible future incremental gains: per-cast integration
with cast averaging (halves drift variance), overlap-constrained integration.

## The "suspiciously quiet surface" question (2026-07-07)

Concern: the inverse's 0–30 m velocities look low. Findings:

1. **Not a solver artifact.** Top-bin values and the drift comparison are insensitive to
   the smoothness weight (0.2/1/5 identical to 3 decimals), bin size (dz 5 vs 10 m), and
   looser surface QC (2-m mask, keep first cell: med|u_top| 0.026 → 0.035, drift
   agreement slightly worse). Near-surface bins do carry ~2× the within-ping noise
   (waves; 0.022 vs 0.009 m/s std) and ~3× fewer samples.
2. **Independent ground truth agrees.** Surface GPS drift vs 0–30 m inverse bins:
   med|Δ| = 0.040 m/s, bias −0.023, across 126 yos. The drift itself is weak: median
   0.038 m/s (p90 0.24).
3. **The surface flow was persistent, not inertial.** Rotary autocorrelation of the
   545 drift vectors: |ρ| = 0.5–0.8 out to 25 h with a weak clockwise drift (−20° at
   the 12.8-h inertial period, not the −360° a rotating inertial current would give).
   November basin interior with a deep mixed layer: the quiet segment-mean surface is
   real at surfacing times.
4. **Honest caveat (the valid core of the intuition):** all our products — inverse,
   shear, DAC, and the drift matching — are yo-segment means. Currents with timescales
   shorter than a yo (storm-forced inertial bursts; drift p90 shows 0.2+ m/s events)
   vector-average toward zero and are under-represented. Resolving them needs a
   time-resolved estimator (per-cast solves are noise-limited in the top bins with only
   ~5–10 pings; a rotary/Kalman formulation à la Stevens-Haas et al. 2022 is the proper
   future tool).

## Task 1 (2026-07-08): shear-content comparison between the methods

New products: [`solve_shear_profile`] (the shear method's pre-integration binned shear)
and [`inverse_shear`] (centered differences of any solver profile — the shear implied by
the inverse). Comparing them on M38 (post-QC, post-calibration, dz = 10 m):

| comparison | r_u | notes |
|---|---|---|
| per-yo bins, direct vs inverse-implied | 0.20 (0.63 at z < 100 m) | rms 1.9×10⁻³ s⁻¹ ≈ the internal-wave shear level |
| same, scale-matched (both profiles 20-m centered-differenced) | 0.22 | scale is not the cause |
| **dive vs climb of the DIRECT product (same method, same yo)** | **0.08** | the reproducibility ceiling |
| mission-mean profile, v | 0.39 (rms 1.6×10⁻⁴; stds 1.4/1.5×10⁻⁴) | agrees |
| mission-mean profile, u | −0.04 (rms 6.7×10⁻⁴) | see below |

**Interpretation.** Per-yo, bin-scale (10–20 m) shear in this environment is dominated
by internal-wave finestructure that decorrelates in the hours between a dive and its
climb — the direct product cannot even reproduce *itself* cast-to-cast (r = 0.08), so
no two estimators can be expected to agree per bin. The expectation "the shear content
should be comparable" holds exactly where signal exceeds the wave noise: the upper
100 m (r = 0.63), the integrated velocity profiles, and the mission-mean v profile.

**The mission-mean u asymmetry** (inverse std 5.3×10⁻⁴ vs direct 1.6×10⁻⁴): the
inverse retains the real sub-inertial shear of the early-mission eddy transit
(magnitudes match the raw-data tilt checks of −5…−6×10⁻⁴ s⁻¹), while the direct
product does not. Working hypothesis: the shear-bias calibration operates entirely
within the ping window, so *track-aligned real mean shear* from the predominantly
eastbound eddy transit leaks into b(r) and is subtracted from the direct product — the
documented heading-diversity caveat, here quantified at O(1–5)×10⁻⁴ s⁻¹ for this
mission's u component — whereas the inverse's shear comes from cross-ping bin structure
and is immune. Task 3 (the 300-m shear feature) will adjudicate the inverse's deep
structure independently via hydrography/thermal wind.

**Guidance** (also in the tutorial): use `inverse_shear` for deterministic/sub-inertial
shear; use `solve_shear_profile` for internal-wave/finestructure statistics (shear
variance) and QC; on low-heading-diversity legs expect the bias calibration to absorb
some real track-aligned mean shear from the *direct* product only.

## Task 3 (2026-07-08): the "strong shear at 300 m" — solved, and a defect found

**Verdict: the feature was an artifact of false bottom-track locks, now screened out by
default.** The investigation, in order:

1. **Localization.** The strongest inverse shear concentrated at 275–415 m, mostly in
   the slope-transit yos (peaks to 1.2×10⁻² s⁻¹ at ~295–395 m), with a milder
   mission-wide elevation in the same band.
2. **Sensitivity.** Smoothness weight: no effect. Referencing: decisive —
   DAC-only median |shear| at 200–450 m in the transit yos is 4.4×10⁻⁴ s⁻¹;
   adding bottom track raises it 10× to 4.2×10⁻³.
3. **Hydrography veto.** The deep pycnocline sits at 160–260 m with
   dσ₀/dz ≤ 0.16 kg/m³ per 100 m (N ≈ 4×10⁻³ s⁻¹): a sustained 10⁻² s⁻¹ shear there
   would put Ri < ¼. Yo-pair thermal wind at 250–400 m supports only
   ~2–5×10⁻⁵ s⁻¹. The measured relative velocities (raw tilt checks) show the smooth
   DAC-only structure, not the sharp feature.
4. **Root cause.** The BT-anchored "absolute" deep velocities were ≈ 0.00 m/s in every
   transit yo while DAC and DAC-only inverse agreed at +0.2…+0.5 — the signature of a
   target moving **with the water**. Classification of all 15,432 three-beam BT locks
   (impossible-bathymetry test: platform later dove deeper than the implied bottom):
   **99.7 % false**, locking on a persistent target 0.6–2.8 m below the transducer
   (median 1.7 m) — near-field/wake, not seafloor (basin depth 1–3 km; the seafloor was
   never within the 30-m BT range on this mission). Anchoring pings to a water-frame
   target contradicts the earth-frame DAC; least squares dumps the contradiction into
   the ocean profile as spurious shear near the top of the anchored depth range — and
   pushes unanchored mid-depth bins to |u| up to 1.6 m/s (the previously flagged
   yos-8–18 outliers, now fully explained; their dive/climb consistency was the blind
   spot — both casts shared the same false anchors).
5. **Fix.** `bt_valid` now defaults to `min_range = 5 m` (genuine seafloor approaches
   have O(10 m) ranges; the false cluster is 0.6–2.8 m) plus an impossible-bathymetry
   screen (implied bottom vs deepest platform record within ±2 h). On M38: 0 of 9,383
   fixes survive (correct — there were no real locks), the transit shear collapses to
   the DAC-only value, the |u| > 0.8 outliers drop from 289 bins to 0, and the
   basin-wide 250–350 m elevation disappears (4.9×10⁻⁴, smoothly decaying with depth —
   **the entire "300-m shear" was this artifact**).

**Amendments to earlier claims (honesty pass):**
- "Independent validation: DAC-only inverse u_g vs unseen bottom track r_v = 0.97,
  med|Δ| ≈ 7 cm/s" — the BT reference was the water-frame target, so this validated the
  **through-water** velocity content of the chain (still a real geometry check; the
  7 cm/s offset ≈ the water speed), *not* an absolute over-ground reference.
- "BT-anchored absolute deep velocities agree to 1.6 cm/s" — self-consistency of
  water-frame anchoring, not an absolute check. Superseded.
- The BT **sign-convention** check (w vs dP/dt, r = 0.90) survives: vertical motion
  relative to a water-borne target still validates the sign chain.
- Headline metrics that never involved BT stand unchanged: dive/climb r = 0.98
  (2 cm/s), DAC closure 5 mm/s, surface drift 4 cm/s, raw-data tilt checks, all
  synthetic and gliderad2cp parity results.
- Task 1's mission-mean u-shear asymmetry: **was BT-injected**, not calibration
  absorption (that hypothesis is retracted) — with hardened screens the DAC-only
  inverse mission-mean u-shear std is 6.8×10⁻⁵ s⁻¹, consistent with the direct
  product.

**Physics answer to the original question:** there is no strong sub-inertial shear
layer at 300 m in this mission. The real deep pycnocline (AW base) sits at 160–260 m
with thermal-wind shear of order 10⁻⁵–10⁻⁴ s⁻¹, and the transit-eddy velocity
structure is smooth across 200–450 m (≈4×10⁻⁴ s⁻¹), exactly as the DAC-only inverse
and the raw relative velocities show.

## Task 4 (2026-07-08): section figures, w products — and a major revision

`examples/m38_currents.jl` now produces four-panel U/V sections (inverse and shear,
both components), a two-panel **w** figure via the new [`solve_w`] product
(`:direct` = binned `U_rel + dP/dt`; `:inverse` = ocean-w bins solved jointly with
per-ping glider w anchored to the pressure-derived vertical velocity — the vertical
analog of a bottom-track constraint), and data-driven symmetric color ranges (99th
percentile of the finite values; M38: ±0.5 m/s horizontal, ±0.04 m/s vertical). The two
w estimates corroborate each other closely; near-surface bins remain wave-contaminated
as documented.

**Revision of the method-intercomparison numbers.** With the false bottom-track locks
screened out (Task 3), the M38 full-mission comparison becomes:

    shear vs inverse:  r_u = 0.977, r_v = 0.979, rms ≈ 0.036 m/s   (was 0.65 / 0.20)
    DAC closure:       median 0.001 m/s                            (was 0.005)

Most of the disagreement previously attributed to shear-method integration noise was in
fact the **inverse being distorted by the false BT anchors**. The corrected picture:
the two methods agree to ≈3–4 cm/s rms mission-wide. What survives of the earlier
narrative: the shear method's per-yo random-walk drift and its ≈+0.09 m/s near-surface
bias against GPS drift (both measured BT-free) are real but second-order; the
cast-to-cast shear-reproducibility ceiling (r = 0.08, method-internal) also stands.
Task-1's per-yo bin-shear comparisons used the contaminated inverse and would improve
somewhat if recomputed, but remain bounded by that sampling ceiling.

## Task 5 (2026-07-08): real-time vs delayed-mode products

Full-pipeline comparison of the two AD2CP data routes, everything else held identical
(gli.sub nav, pld1.sub CTD sound speed, QC thresholds, IGRF declination, shear-bias
calibration, DAC): the **real-time** `$PNOR` ASCII telemetry stream (`load_pnor`;
0.01 m/s velocity quantization, 0.1° attitude, no accelerometer, no BT records) vs the
**delayed-mode** full-resolution `.ad2cp` binary. Script:
`examples/m38_realtime_vs_delayed.jl`; gated acceptance test in the suite.

Coverage first (the Task-6 machinery reports it directly): the stream carries
123,950 of 124,752 ensembles (99.4%). The payload stopped writing the stream on
2022-11-27 — segment sequence complete, `missing_segments` empty, so this is the
payload configuration, not transfer loss — while the instrument kept recording sparse
internal bursts through 2023-03-01 (750 ensembles, three months, delayed-only).
The real-time record covers the entire main pinging period.

Product agreement over the 127 yos both sides solve (identical yo sets — real-time
loses no segments; common (yo, z) bins, nobs > 10):

    inverse u:  r = 0.9996   rms = 4.6 mm/s   bias  0.0 mm/s
    inverse v:  r = 0.9997   rms = 4.1 mm/s   bias  0.0 mm/s
    shear u:    r = 0.9897   rms = 24.9 mm/s  bias +0.1 mm/s
    shear v:    r = 0.9911   rms = 22.9 mm/s  bias -0.1 mm/s
    w (direct): r = 0.9603   rms = 2.6 mm/s   bias  0.0 mm/s

Depth structure: the inverse difference is flat at 4–6 mm/s from the surface to
1000 m (5.7 mm/s in 600–1000 m, where fewer samples average the quantization noise);
the shear-method difference runs 18–31 mm/s because vertical integration accumulates
the quantized-sample noise that the inverse localizes per bin.

**Conclusions.** (1) A real-time/onboard product built from the telemetry stream is
essentially the delayed product: the inverse-method penalty is ~5 mm/s rms with zero
bias — an order of magnitude below the ~3–4 cm/s method/sampling uncertainty
established in Task 3. The 0.01 m/s per-sample quantization averages down as expected
(~200 samples per bin → few mm/s). (2) The stream's missing pieces (accelerometer →
pass `look` explicitly; BT records → irrelevant on M38, which has no genuine BT;
magnetometer → declination comes from nav anyway) cost nothing here. (3) Prefer the
shear method last in real-time settings: it is the one product measurably degraded
(≈2.5 cm/s rms) by stream quantization.
