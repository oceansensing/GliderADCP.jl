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
