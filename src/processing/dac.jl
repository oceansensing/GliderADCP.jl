# Layer 3 — depth-averaged current (DAC) and surface drift from navigation.
#
# SeaExplorer principle: while submerged (`DeadReckoning == 1`) the vehicle dead-reckons
# its position from heading + a through-water flight model, ignoring currents. At
# surfacing, the first GPS fix (`DeadReckoning == 0`) snaps the position; the jump
# between the last dead-reckoned position and that fix is the current-induced
# displacement accumulated over the submerged interval:
#
#     DAC = (pos_fix − pos_DR_end) / (t_fix − t_submerged_start)
#
# This is the standard glider DAC (Rudnick et al. 2018; Gradone et al. 2023 report
# 1–2 cm/s RMS accuracy for the equivalent Slocum estimate).

const _EARTH_R = 6.371e6  # m

"""
    lonlat_to_dxdy(lon0, lat0, lon1, lat1) -> (dx, dy)

Local-tangent displacement in meters (east, north) from position 0 to position 1
(decimal degrees; spherical earth, cosine of the mean latitude).
"""
function lonlat_to_dxdy(lon0, lat0, lon1, lat1)
    latm = (lat0 + lat1) / 2
    dx = deg2rad(lon1 - lon0) * _EARTH_R * cosd(latm)
    dy = deg2rad(lat1 - lat0) * _EARTH_R
    return dx, dy
end

"""
    compute_dac(nav::GliderNav; min_duration=600.0, min_depth=10.0,
                max_speed=1.5, max_fix_delay=900.0) -> DataFrame

Depth-averaged current per submerged segment from SeaExplorer navigation.

Each maximal `DeadReckoning == 1` block bounded by finite-position GPS fixes yields one
estimate. Quality control drops segments that are too short (`min_duration`, s), too
shallow (`min_depth`, m — surface drift intervals), whose first fix arrives too long
after the last DR record (`max_fix_delay`, s), or that imply unphysical speeds
(`max_speed`, m/s).

Returns a `DataFrame` with one row per accepted segment:
`yo, t_start, t_end, t_mid, duration, lon0, lat0, lon_dr, lat_dr, lon_fix, lat_fix,
maxdepth, u, v` — where `(u, v)` is the DAC (m/s, east/north), `(lon0, lat0)` the last
fix before diving, `(lon_dr, lat_dr)` the final dead-reckoned position and
`(lon_fix, lat_fix)` the first fix after surfacing.
"""
function compute_dac(nav::GliderNav; min_duration::Real=600.0, min_depth::Real=10.0,
                     max_speed::Real=1.5, max_fix_delay::Real=900.0)
    n = length(nav)
    dr = nav.deadreckoning
    isfix(i) = dr[i] == 0 && isfinite(nav.lon[i]) && isfinite(nav.lat[i])

    rows = NamedTuple[]
    yo = 0
    i = 1
    while i <= n
        if dr[i] != 1
            i += 1
            continue
        end
        # maximal submerged block [i1, i2]
        i1 = i
        i2 = i
        while i2 < n && dr[i2+1] == 1
            i2 += 1
        end
        i = i2 + 1

        # last fix before, first fix after
        ib = 0
        for j in i1-1:-1:1
            isfix(j) && (ib = j; break)
        end
        ia = 0
        for j in i2+1:n
            isfix(j) && (ia = j; break)
        end
        (ib == 0 || ia == 0) && continue
        isfinite(nav.lon[i2]) && isfinite(nav.lat[i2]) || continue

        # DR error (and current drift) accrues from the last pre-dive fix to the first
        # post-surfacing fix — the displacement window is fix-to-fix
        duration = nav.t[ia] - nav.t[ib]
        fix_delay = nav.t[ia] - nav.t[i2]
        depths = filter(isfinite, nav.depth[i1:i2])
        maxdepth = isempty(depths) ? NaN : maximum(depths)

        dx, dy = lonlat_to_dxdy(nav.lon[i2], nav.lat[i2], nav.lon[ia], nav.lat[ia])
        u = dx / duration
        v = dy / duration

        yo += 1
        ok = duration >= min_duration && fix_delay <= max_fix_delay &&
             isfinite(maxdepth) && maxdepth >= min_depth &&
             isfinite(u) && isfinite(v) && abs(u) <= max_speed && abs(v) <= max_speed
        ok || continue

        push!(rows, (
            yo = yo,
            t_start = nav.time[ib], t_end = nav.time[ia],
            t_mid = nav.time[ib] + Millisecond(round(Int, 500duration)),
            duration = duration,
            lon0 = nav.lon[ib], lat0 = nav.lat[ib],
            lon_dr = nav.lon[i2], lat_dr = nav.lat[i2],
            lon_fix = nav.lon[ia], lat_fix = nav.lat[ia],
            maxdepth = maxdepth,
            u = u, v = v,
        ))
    end
    return DataFrame(rows)
end

"""
    surface_drift(nav::GliderNav; min_gap=30.0, max_gap=1200.0, max_speed=1.5) -> DataFrame

Near-surface drift velocities from consecutive GPS fixes within the same surface
interval (no submerged records between them). Returns `t_mid, duration, lon, lat, u, v`
per fix pair — a near-surface velocity constraint for the inverse solution.
"""
function surface_drift(nav::GliderNav; min_gap::Real=30.0, max_gap::Real=1200.0,
                       max_speed::Real=1.5)
    n = length(nav)
    rows = NamedTuple[]
    prev = 0
    for i in 1:n
        if nav.deadreckoning[i] == 1
            prev = 0
        elseif nav.deadreckoning[i] == 0 && isfinite(nav.lon[i]) && isfinite(nav.lat[i])
            if prev > 0
                dt = nav.t[i] - nav.t[prev]
                if min_gap <= dt <= max_gap
                    dx, dy = lonlat_to_dxdy(nav.lon[prev], nav.lat[prev], nav.lon[i], nav.lat[i])
                    u, v = dx / dt, dy / dt
                    if abs(u) <= max_speed && abs(v) <= max_speed
                        push!(rows, (
                            t_mid = nav.time[prev] + Millisecond(round(Int, 500dt)),
                            duration = dt,
                            lon = (nav.lon[prev] + nav.lon[i]) / 2,
                            lat = (nav.lat[prev] + nav.lat[i]) / 2,
                            u = u, v = v,
                        ))
                    end
                end
            end
            prev = i
        end
    end
    return DataFrame(rows)
end
