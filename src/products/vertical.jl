# Layer 5 — vertical water velocity (flight-model-free).
#
# With an ADCP the vertical water velocity comes directly from the measured relative
# vertical velocity plus the glider's own vertical motion from the pressure record:
#     w_water(z, t) = U_rel(z, t) + w_glider(t),   w_glider = −d(depth)/dt  (+up)
# No Merckelbach-style flight model is required for w when an ADCP is present (the
# model remains useful for angle-of-attack/performance studies; cf. Merckelbach et al.
# 2010/2019, Todd et al. 2017 §flight constraints).

"""
    vertical_velocity(p::ProcessedPings; max_gap=60.0) -> Matrix

Absolute vertical water velocity (m/s, positive up) on the ping × offset grid:
`w = U_rel + w_glider`, with `w_glider = −d(depth)/dt` by centered differences
(ping gaps larger than `max_gap` seconds yield NaN).
"""
function vertical_velocity(p::ProcessedPings; max_gap::Real=60.0)
    n = length(p)
    wg = fill(NaN, n)
    for i in 2:n-1
        dt = p.t[i+1] - p.t[i-1]
        (0 < dt <= 2max_gap && isfinite(p.depth[i+1]) && isfinite(p.depth[i-1])) || continue
        wg[i] = -(p.depth[i+1] - p.depth[i-1]) / dt
    end
    return p.U .+ wg'
end
