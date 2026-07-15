# Layer 1 — Slocum glider data ingestion (platform-agnostic tables).
#
# Works from any Slocum-derived table (ERDDAP tabledap export, Python dbdreader, or the
# user's pure-Julia SlocumIO.jl — https://github.com/oceansensing/SlocumIO.jl):
#
#   using SlocumIO
#   d = MultiDBD(dir="…", eng=true, sci=true)
#   df = DataFrame(get_sync(d, "m_water_vx", "m_water_vy", "m_gps_mag_var",
#                              "m_lat", "m_lon", "m_depth", ...))
#
# then `slocum_nav(df)` / `dac_from_slocum(df)` feed the standard pipeline.

_col(df, names...) = begin
    for n in names
        hasproperty(df, n) && return df[!, n]
    end
    nothing
end

# Per-row DateTime, tolerant of `missing` and NaN stamps (ERDDAP exports routinely
# carry a few): invalid rows map to the 1970 epoch sentinel and are dropped by the
# callers with a warning — never an InexactError/MethodError crash.
const _SLOCUM_EPOCH = DateTime(1970)
_to_datetime(t::AbstractVector{<:DateTime}) = collect(t)
_to_datetime(t::AbstractVector) =
    DateTime[x isa DateTime ? x :
             (x !== missing && isfinite(Float64(x))) ? unix2datetime(Float64(x)) :
             _SLOCUM_EPOCH for x in t]
_valid_time(x::DateTime) = x > DateTime(1971)

# valid numeric cell under BOTH conventions: `missing` (ERDDAP/CSV) and NaN
# (dbdreader / SlocumIO gap fills)
_oknum(x) = x !== missing && isfinite(x)

"""
    slocum_nav(df) -> GliderNav

Build a [`GliderNav`](@ref) from a Slocum table. Recognized columns (first match wins):
time (`time`), position (`latitude`/`m_gps_lat`, `longitude`/`m_gps_lon`, decimal deg),
depth (`depth`/`m_depth`), attitude (`m_heading`/`m_pitch`/`m_roll`, radians),
declination (`m_gps_mag_var`, radians). Missing columns become NaN;
`deadreckoning`/`navstate` are set to unknown (Slocum DAC comes from
[`dac_from_slocum`](@ref) instead of DR/GPS jumps).
"""
function slocum_nav(df::DataFrame)
    tcol = _col(df, :time)
    tcol === nothing && error("slocum_nav: no `time` column")
    tall = _to_datetime(tcol)
    keep = findall(_valid_time, tall)
    ndrop = nrow(df) - length(keep)
    ndrop > 0 && @warn "slocum_nav: dropped $ndrop row(s) with missing/epoch time"
    dfk = df[keep, :]
    time = tall[keep]
    n = length(time)
    f(names...; scale=1.0) = begin
        c = _col(dfk, names...)
        c === nothing ? fill(NaN, n) :
            Float64[x === missing ? NaN : Float64(x) for x in c] .* scale
    end
    GliderNav(time, datetime2unix.(time),
        f(:longitude, :m_gps_lon, :lon), f(:latitude, :m_gps_lat, :lat),
        f(:m_heading, :heading; scale=180 / π),
        f(:m_gps_mag_var; scale=180 / π),
        f(:m_pitch, :pitch; scale=180 / π), f(:m_roll, :roll; scale=180 / π),
        f(:depth, :m_depth), fill(Int16(-1), n), fill(Int8(-1), n), fill(NaN, n), dfk)
end

"""
    dac_from_slocum(df; by=:source_file, min_depth=10.0) -> DataFrame

Per-segment depth-averaged current from the glider's own dead-reckoned estimate
(`m_water_vx/vy`, magnetic frame), rotated to true east/north by `m_gps_mag_var`
(Gradone et al. 2023 recipe: last non-missing value per segment, mean declination).
Segments come from the `by` column (Slocum `source_file` ≈ surfacing-to-surfacing).
Output matches the [`compute_dac`](@ref) schema used by the solvers
(`yo, t_start, t_end, t_mid, duration, u, v`).
"""
function dac_from_slocum(df::DataFrame; by::Symbol=:source_file, min_depth::Real=10.0)
    hasproperty(df, by) || error("dac_from_slocum: no `$by` column")
    tcol = _col(df, :time)
    tcol === nothing && error("dac_from_slocum: no `time` column")
    time = _to_datetime(tcol)
    rows = NamedTuple[]
    yo = 0
    for g in groupby(DataFrame(df; copycols=false), by)
        idx = parentindices(g)[1]
        t = filter(_valid_time, time[idx])
        isempty(t) && continue
        vx = _col(g, :m_water_vx); vy = _col(g, :m_water_vy)
        (vx === nothing || vy === nothing) && continue
        iv = findlast(i -> _oknum(vx[i]) && _oknum(vy[i]), 1:nrow(g))
        iv === nothing && continue
        dep = _col(g, :depth, :m_depth)
        if dep !== nothing
            dmax = maximum(Float64[x for x in dep if _oknum(x)]; init=0.0)
            dmax < min_depth && continue
        end
        mv = _col(g, :m_gps_mag_var)
        mvok = mv === nothing ? Float64[] : Float64[x for x in mv if _oknum(x)]
        mvdeg = isempty(mvok) ? 0.0 : rad2deg(mean(mvok))
        u0, v0 = Float64(vx[iv]), Float64(vy[iv])
        c, s = cosd(mvdeg), sind(mvdeg)
        yo += 1
        t1, t2 = extrema(t)
        push!(rows, (yo=yo, t_start=t1, t_end=t2,
            t_mid=t1 + Millisecond(round(Int, (t2 - t1).value / 2)),
            duration=(t2 - t1).value / 1000,
            u=u0 * c - v0 * s, v=u0 * s + v0 * c))
    end
    isempty(rows) && return DataFrame(yo=Int[], t_start=DateTime[], t_end=DateTime[],
        t_mid=DateTime[], duration=Float64[], u=Float64[], v=Float64[])
    return sort!(DataFrame(rows), :t_start)
end
