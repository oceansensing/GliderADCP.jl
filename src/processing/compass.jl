# Layer 2 diagnostics — compass/magnetometer health.
#
# The magnitude of the magnetometer vector should be independent of heading and pitch
# for a clean, well-calibrated compass; heading-dependent |B| indicates hard/soft-iron
# contamination that rotates ENU velocities (von Appen 2015). This is a *diagnostic*
# (detection and flagging); a deviation correction is a research task.

"""
    compass_field_check(adcp; nsector=12) -> (table, ptp_fraction)

Median total magnetometer field |B| per heading sector (counts), overall
peak-to-peak variation as a fraction of the mission median, and sample counts.
`ptp_fraction` ≳ 0.05 suggests heading-dependent field contamination worth
investigating before trusting absolute current directions to a few degrees.
"""
function compass_field_check(a::AD2CPData; nsector::Int=12)
    Bm = vec(sqrt.(sum(abs2, a.mag; dims=1)))
    sect = zeros(Int, length(a))
    for i in 1:length(a)
        h = a.heading[i]
        sect[i] = isfinite(h) && isfinite(Bm[i]) ?
                  clamp(floor(Int, mod(h, 360) / (360 / nsector)) + 1, 1, nsector) : 0
    end
    rows = NamedTuple[]
    meds = Float64[]
    for s in 1:nsector
        idx = findall(==(s), sect)
        length(idx) < 50 && continue
        m = median(Bm[idx])
        push!(meds, m)
        push!(rows, (sector=s, heading=(s - 0.5) * 360 / nsector, medB=m, n=length(idx)))
    end
    tbl = DataFrame(rows)
    ptp = isempty(meds) ? NaN : (maximum(meds) - minimum(meds)) / median(meds)
    return tbl, ptp
end
