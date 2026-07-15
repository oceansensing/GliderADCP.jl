# Layer 2/3 calibration — range-dependent ("shear") bias.
#
# Glider AD2CPs exhibit a small systematic decay of measured relative velocity with
# range from the transducer (±few mm/s across the profiling window; Todd et al. 2017;
# gliderad2cp `process_bias`). It is invisible to per-sample QC, but the shear method
# integrates it into a spurious profile tilt of O(0.1–0.2 m/s per 500 m). On the M38
# reference mission the measured along-track slope is −3.2×10⁻⁴ s⁻¹, identical on dives
# and climbs, with a null cross-track component (docs/research/m38_validation.md).
#
# Calibration: the mission-mean *within-ping* anomaly of relative velocity vs offset,
# in the glider track frame. Real ocean shear projects onto the track frame with the
# heading and averages out over a mission with varied headings; the instrument bias is
# track-locked and survives. Subtracting the calibrated profile is mean-free per ping,
# so it changes vertical structure only (DAC/glider-velocity content is untouched).

# through-water speed proxy for a ping: magnitude of its mean relative velocity
function _ping_speed(p::ProcessedPings, i::Int, fin::Vector{Int})
    se = 0.0; sn = 0.0
    for k in fin
        se += p.E[k, i]; sn += p.N[k, i]
    end
    return hypot(se, sn) / length(fin)
end

"""
    shear_bias(p::ProcessedPings; min_count=1000, velocity_scaled=false) -> NamedTuple

Mission-calibrated range-dependent bias profile in the glider track frame. Returns
`(offsets, along, cross, n, slope_along, slope_cross, heading_concentration,
velocity_scaled, ref_speed)`.

The calibration averages the *adjacent-pair differences* of the track-frame relative
velocities — the exact sample population the shear estimator consumes — and integrates
them into mean-removed bias profiles `along`/`cross` (m/s on `p.offsets`; pairs with
fewer than `min_count` samples contribute zero difference). With `velocity_scaled=true`
the bias is modeled proportional to each ping's through-water speed (the
speed-correlated mechanism of Todd et al. 2017 / gliderad2cp's velocity-dependent
correction): `along`/`cross` hold dimensionless coefficient profiles `c(k)` such that
the bias of ping *i* is `c(k) · speedᵢ`; `slope_*` are reported in s⁻¹ at the mission
median speed `ref_speed` for comparability.

`heading_concentration` is the mission heading resultant `R ∈ [0, 1]`; values near 1
mean a one-heading mission where real ocean shear can leak into the calibration
(warned above 0.8).
"""
function shear_bias(p::ProcessedPings; min_count::Int=1000, velocity_scaled::Bool=false)
    nk = length(p.offsets)
    sa = zeros(nk); sc = zeros(nk)
    wa = zeros(nk)                    # Σ speed² for the scaled fit, Σ 1 otherwise
    n = zeros(Int, nk)
    shs = 0.0; chs = 0.0; nh = 0
    speeds = Float64[]
    at = Float64[]; ct = Float64[]; fin = Int[]
    for i in 1:length(p)
        h = p.heading[i]
        isfinite(h) || continue
        sh, ch = sind(h), cosd(h)
        empty!(fin)
        for k in 1:nk
            (isfinite(p.E[k, i]) && isfinite(p.N[k, i])) && push!(fin, k)
        end
        length(fin) < 5 && continue
        spd = velocity_scaled ? _ping_speed(p, i, fin) : 1.0
        (isfinite(spd) && spd > 0) || continue
        velocity_scaled && push!(speeds, spd)
        # calibrate on ADJACENT-PAIR differences — the exact population the shear
        # estimator consumes (with partial coverage, mean-of-differences ≠
        # difference-of-means, so a per-offset mean profile under-corrects)
        for k in 1:nk-1
            (isfinite(p.E[k, i]) && isfinite(p.N[k, i]) &&
             isfinite(p.E[k+1, i]) && isfinite(p.N[k+1, i])) || continue
            dat = (p.E[k+1, i] - p.E[k, i]) * sh + (p.N[k+1, i] - p.N[k, i]) * ch
            dct = -(p.E[k+1, i] - p.E[k, i]) * ch + (p.N[k+1, i] - p.N[k, i]) * sh
            # LS for diff = δ(k)·spd:  δ(k) = Σ(diff·spd)/Σ(spd²);  spd ≡ 1 unscaled
            sa[k] += dat * spd
            sc[k] += dct * spd
            wa[k] += spd^2
            n[k] += 1
        end
        shs += sh; chs += ch; nh += 1
    end
    # pair-difference biases → integrated (mean-removed) bias profiles B(k)
    δa = [n[k] >= min_count ? sa[k] / wa[k] : 0.0 for k in 1:nk-1]
    δc = [n[k] >= min_count ? sc[k] / wa[k] : 0.0 for k in 1:nk-1]
    along = vcat(0.0, cumsum(δa))
    cross = vcat(0.0, cumsum(δc))
    along .-= mean(along)
    cross .-= mean(cross)
    if velocity_scaled && isempty(speeds)
        @warn "shear_bias: velocity_scaled requested but no ping qualified for a " *
              "speed estimate — returning the unscaled estimate"
    end
    ref_speed = velocity_scaled && !isempty(speeds) ? median(speeds) : 1.0
    fitslope(b) = begin
        g = findall(isfinite, b)
        length(g) > 3 ? ref_speed * cov(p.offsets[g], b[g]) / var(p.offsets[g]) : NaN
    end
    R = nh > 0 ? hypot(shs, chs) / nh : NaN
    isfinite(R) && R > 0.8 &&
        @warn "shear_bias: headings are concentrated (R=$(round(R, digits=2))) — real " *
              "ocean shear may leak into the calibration; interpret with care"
    return (offsets=copy(p.offsets), along=along, cross=cross, n=n,
        slope_along=fitslope(along), slope_cross=fitslope(cross),
        heading_concentration=R, velocity_scaled=velocity_scaled, ref_speed=ref_speed)
end

"""
    apply_shear_bias!(p::ProcessedPings, b; components=:both) -> p

Subtract the calibrated track-frame bias profile from the per-ping relative velocities
(rotated to ENU with each ping's heading). The subtracted profile is re-demeaned over
each ping's finite cells, so **ping-mean velocities are unchanged by construction** —
the correction alters vertical structure only (the inverse's glider-velocity and DAC
content is untouched). `components = :along` limits the correction to the along-track
profile (the cross-track bias is typically negligible).
"""
function apply_shear_bias!(p::ProcessedPings, b; components::Symbol=:both)
    nk = length(p.offsets)
    length(b.along) == nk || error("apply_shear_bias!: offset grids differ")
    ba = [isfinite(x) ? x : 0.0 for x in b.along]
    bc = components === :both ? [isfinite(x) ? x : 0.0 for x in b.cross] : zeros(nk)
    scaled = get(b, :velocity_scaled, false)
    fin = Int[]
    for i in 1:length(p)
        h = p.heading[i]
        isfinite(h) || continue
        sh, ch = sind(h), cosd(h)
        empty!(fin)
        for k in 1:nk
            (isfinite(p.E[k, i]) && isfinite(p.N[k, i])) && push!(fin, k)
        end
        isempty(fin) && continue
        spd = scaled ? _ping_speed(p, i, fin) : 1.0
        (isfinite(spd) && spd > 0) || continue
        mba = mean(ba[k] for k in fin)
        mbc = mean(bc[k] for k in fin)
        for k in fin
            da = (ba[k] - mba) * spd
            dc = (bc[k] - mbc) * spd
            p.E[k, i] -= da * sh - dc * ch
            p.N[k, i] -= da * ch + dc * sh
        end
    end
    return p
end

"""
    calibrate_shear_bias!(p::ProcessedPings; passes=3, kwargs...) -> Vector

Iterated estimate-and-subtract of the range-dependent bias (partial-coverage pings make
a single pass incomplete). Returns the fitted along-track slope (s⁻¹) after each pass —
the last value is the residual. `kwargs` are forwarded to [`shear_bias`](@ref)
(e.g. `velocity_scaled=true`).
"""
function calibrate_shear_bias!(p::ProcessedPings; passes::Int=3, kwargs...)
    slopes = Float64[]
    for _ in 1:passes
        b = shear_bias(p; kwargs...)
        push!(slopes, b.slope_along)
        apply_shear_bias!(p, b)
    end
    b = shear_bias(p; kwargs...)
    push!(slopes, b.slope_along)
    return slopes
end
