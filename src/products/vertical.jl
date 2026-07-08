# Layer 5 — vertical water velocity (flight-model-free).
#
# With an ADCP the vertical water velocity comes directly from the measured relative
# vertical velocity plus the glider's own vertical motion from the pressure record:
#     w_water(z, t) = U_rel(z, t) + w_glider(t),   w_glider = −d(depth)/dt  (+up)
# No Merckelbach-style flight model is required for w when an ADCP is present (the
# model remains useful for angle-of-attack/performance studies; cf. Merckelbach et al.
# 2010/2019, Todd et al. 2017 §flight constraints).

"""
    glider_w(p::ProcessedPings; max_gap=60.0) -> Vector

Glider vertical velocity (m/s, positive up) from the pressure record:
`−d(depth)/dt` by centered differences (gaps larger than `max_gap` s yield NaN).
"""
function glider_w(p::ProcessedPings; max_gap::Real=60.0)
    n = length(p)
    wg = fill(NaN, n)
    for i in 2:n-1
        dt = p.t[i+1] - p.t[i-1]
        (0 < dt <= 2max_gap && isfinite(p.depth[i+1]) && isfinite(p.depth[i-1])) || continue
        wg[i] = -(p.depth[i+1] - p.depth[i-1]) / dt
    end
    return wg
end

"""
    vertical_velocity(p::ProcessedPings; max_gap=60.0) -> Matrix

Absolute vertical water velocity (m/s, positive up) on the ping × offset grid:
`w = U_rel + w_glider`, with `w_glider = −d(depth)/dt` by centered differences
(ping gaps larger than `max_gap` seconds yield NaN).
"""
vertical_velocity(p::ProcessedPings; max_gap::Real=60.0) = p.U .+ glider_w(p; max_gap)'

"""
    solve_w(pings::ProcessedPings, segments::DataFrame;
            method=:direct, dz=10.0, min_bin_obs=4, min_pings=30,
            wanchor=5.0, wsmooth=1.0) -> DataFrame

Vertical-water-velocity profiles per segment, by either of two methods over the same
samples (`segments` needs `yo, t_start, t_end, t_mid`; a [`compute_dac`](@ref) table
works):

  * `:direct` — bin the per-sample absolute w (= `U_rel + w_glider`,
    [`vertical_velocity`](@ref)) by cell depth: the measurement route.
  * `:inverse` — the inverse machinery applied to the vertical component: ocean-w depth
    bins and per-ping glider w solved jointly, with the glider-w unknowns anchored
    (weight `wanchor`) to the pressure-derived `−d(depth)/dt` — the vertical analog of
    a bottom-track constraint — plus `wsmooth` second-difference smoothing.

Returns `yo, t_mid, z, w, nobs` (w in m/s, positive up).
"""
function solve_w(p::ProcessedPings, segments::DataFrame;
                 method::Symbol=:direct, dz::Real=10.0, min_bin_obs::Int=4,
                 min_pings::Int=30, wanchor::Real=5.0, wsmooth::Real=1.0)
    out = DataFrame(yo=Int[], t_mid=DateTime[], z=Float64[], w=Float64[], nobs=Int[])
    wg = glider_w(p)
    wabs = method === :direct ? p.U .+ wg' : nothing
    for row in eachrow(segments)
        idx = segment_indices(p, row.t_start, row.t_end)
        length(idx) >= min_pings || continue
        if method === :direct
            acc = Dict{Int,Vector{Float64}}()
            for i in idx, k in 1:length(p.offsets)
                (isfinite(wabs[k, i]) && isfinite(p.celldepth[k, i]) &&
                 p.celldepth[k, i] >= 0) || continue
                push!(get!(acc, floor(Int, p.celldepth[k, i] / dz), Float64[]), wabs[k, i])
            end
            for kb in sort(collect(keys(acc)))
                v = acc[kb]
                length(v) >= min_bin_obs || continue
                push!(out, (row.yo, row.t_mid, (kb + 0.5) * dz, mean(v), length(v)))
            end
        elseif method === :inverse
            fin = [i for i in idx if isfinite(wg[i])]
            length(fin) >= min_pings || continue
            anchor = DataFrame(t=p.t[fin], u=wg[fin], v=zeros(length(fin)))
            gd = filter(isfinite, p.depth[idx])
            isempty(gd) && continue
            sol = invert_segment(view(p.U, :, idx), zeros(length(p.offsets), length(idx)),
                view(p.celldepth, :, idx), p.t[idx], maximum(gd);
                bt=anchor, opts=InverseOptions(dz=Float64(dz), wdac=0.0, wbt=Float64(wanchor),
                    wsmooth_ocean=Float64(wsmooth), min_pings=min_pings,
                    min_bin_obs=min_bin_obs))
            sol === nothing && continue
            for k in eachindex(sol.z)
                push!(out, (row.yo, row.t_mid, sol.z[k], sol.u[k], sol.nobs[k]))
            end
        else
            error("solve_w: method must be :direct or :inverse")
        end
    end
    return out
end
