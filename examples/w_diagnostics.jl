# Vertical-velocity quality: the two systematic artifacts, and what the
# calibrations do about them — per mission, both resolved against cell range.
#
#   Artifact A  dive/climb asymmetry growing with range = the vertical projection
#               of the range-dependent beam bias. Shown BEFORE and AFTER
#               calibrate_vertical_bias!; the ocean is symmetric between dives and
#               climbs, so whatever survives that difference is instrumental.
#   Artifact B  a downward bias below |pitch| ≈ 6° = unsteady flight through the
#               inflection. Masked by default (w_min on |glider_w|); the panel
#               shows the raw pitch dependence the screen removes.
#
# Evidence and scope: QA/QC guide §3c, validation doc 2026-07-16 w entry.
#
# Run (CairoMakie comes from the user's @ocean env via the stacked load path;
# the package itself stays plot-free):
#
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/w_diagnostics.jl
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/w_diagnostics.jl m38

using GliderADCP
using CairoMakie
using DataFrames, Dates, Statistics, NaNStatistics
using Printf
include("missions.jl")

const OUT = joinpath(@__DIR__, "output")
mkpath(OUT)

q(x, p) = (v = sort(filter(isfinite, x)); isempty(v) ? NaN :
           v[clamp(round(Int, p * length(v)), 1, length(v))])

const OFFB = [(2, 6), (6, 10), (10, 14), (14, 18), (18, 24), (24, 30)]
const PB = [(0, 3), (3, 6), (6, 10), (10, 14), (14, 18), (18, 22), (22, 26), (26, 40)]

# median(climb) − median(dive) of w per offset band, steady flight only: zero for
# a clean instrument, since the ocean does not care which way the glider points
function dive_climb_asymmetry(p, pitch)
    wc = vertical_velocity(p; w_min=0.0)          # unscreened: measure the raw effect
    out = Float64[]
    for (lo, hi) in OFFB
        ks = findall(o -> lo <= abs(o) < hi, p.offsets)
        if isempty(ks)
            push!(out, NaN); continue
        end
        dv = Float64[]; cv = Float64[]
        for i in 1:length(p)
            (isfinite(p.depth[i]) && p.depth[i] > 20 &&
             isfinite(pitch[i]) && abs(pitch[i]) >= 14) || continue
            v = [wc[k, i] for k in ks if isfinite(wc[k, i])]
            length(v) >= 2 || continue
            pitch[i] < 0 ? push!(dv, mean(v)) : push!(cv, mean(v))
        end
        push!(out, (length(dv) < 50 || length(cv) < 50) ? NaN : q(cv, 0.5) - q(dv, 0.5))
    end
    return out
end

res = Dict{String,Any}()
keys_ = selected_missions()
for key in keys_
    m = MISSIONS[key]
    @info "════════ $(m.label) ════════"
    nav = load_seaexplorer_nav([joinpath(m.dir, "delayed/nav/logs"), joinpath(m.dir, "glimpse")];
                               stream="$(m.prefix).gli.sub")
    lat = round(nanmedian(nav.lat), digits=1)
    adcp = read_ad2cp(joinpath(m.dir, m.binary))
    qc!(adcp)
    p = process_pings(adcp; lat, declination=magnetic_declination(nav, adcp.t))
    calibrate_shear_bias!(p)
    pitch = Float64.(adcp.pitch)                  # aligned 1:1 with the pings

    before = dive_climb_asymmetry(p, pitch)

    # Artifact B: raw pitch dependence (screen off, so the bias is visible)
    wc = vertical_velocity(p; w_min=0.0)
    selmid = findall(o -> 6.0 <= abs(o) <= 20.0, p.offsets)
    wping = fill(NaN, length(p))
    for i in 1:length(p)
        v = [wc[k, i] for k in selmid if isfinite(wc[k, i])]
        length(v) >= 3 && (wping[i] = mean(v))
    end
    deep = [i for i in 1:length(p) if isfinite(p.depth[i]) && p.depth[i] > 20 &&
            isfinite(pitch[i]) && isfinite(wping[i])]
    px = Float64[]; pw = Float64[]
    for (lo, hi) in PB
        ii = [i for i in deep if lo <= abs(pitch[i]) < hi]
        push!(px, (lo + hi) / 2)
        push!(pw, length(ii) < 20 ? NaN : q(wping[ii], 0.5))
    end

    vslopes = calibrate_vertical_bias!(p)
    after = dive_climb_asymmetry(p, pitch)
    @printf("  vertical-bias slope %.3e → residual %.1e s⁻¹\n", vslopes[1], vslopes[end])
    @printf("  asymmetry (mm/s) by offset  before: %s\n",
            join([@sprintf("%+6.1f", 1000x) for x in before], " "))
    @printf("                               after: %s\n",
            join([@sprintf("%+6.1f", 1000x) for x in after], " "))

    dac = compute_dac(nav, p; fallback=flight_model(nav))
    wd = solve_w(p, dac)                          # screened + calibrated product
    wi = solve_w(p, dac; method=:inverse)
    j = innerjoin(wd, wi; on=[:yo, :z], makeunique=true)
    g = (j.nobs .> 10) .& isfinite.(j.w) .& isfinite.(j.w_1)
    @printf("  :direct vs :inverse  n=%d  r=%.4f  rms=%.4f m/s | median w %+.5f m/s\n",
            count(g), cor(j.w[g], j.w_1[g]),
            sqrt(mean((j.w[g] .- j.w_1[g]) .^ 2)), nanmedian(wd.w))
    res[m.label] = (; before, after, px, pw, sec=grid_profiles(wd; fields=(:w, :w)))
end

labels = [MISSIONS[k].label for k in keys_]
cols = [:steelblue, :firebrick, :seagreen, :darkorange]
offx = [(lo + hi) / 2 for (lo, hi) in OFFB]

fig = Figure(size=(1500, 900))
Label(fig[0, 1:2], "Vertical velocity: the two systematic artifacts and their corrections";
      fontsize=21, font=:bold)

sec = res[labels[1]].sec
crw = ceil(q(vec(sec.U), 0.99) * 200) / 200
ax0 = Axis(fig[1, 1:2]; ylabel="depth (m)", yreversed=true,
    title="$(labels[1]) w section — calibrated and inflection-screened product",
    xticks=(let n = length(sec.t), xt = round.(Int, range(1, n; length=8))
                (xt, Dates.format.(sec.t[xt], "dd u")) end))
hm = heatmap!(ax0, 1:length(sec.t), sec.z, permutedims(sec.U);
    colormap=:balance, colorrange=(-crw, crw))
Colorbar(fig[1, 3], hm; label="w (m/s)")

ax1 = Axis(fig[2, 1]; xlabel="cell offset from glider (m)",
    ylabel="climb − dive median w (m/s)",
    title="Artifact A: dive/climb asymmetry vs range\nopen = before, filled = after calibrate_vertical_bias!")
for (i, L) in enumerate(labels)
    r = res[L]
    scatterlines!(ax1, offx, r.before; color=(cols[i], 0.45), markersize=10,
        marker=:circle, linestyle=:dash, label="$L before")
    scatterlines!(ax1, offx, r.after; color=cols[i], markersize=11, label="$L after")
end
hlines!(ax1, [0.0]; color=:gray, linestyle=:dash)
axislegend(ax1; position=:lb, nbanks=2, labelsize=9)

ax2 = Axis(fig[2, 2]; xlabel="|pitch| (°)", ylabel="median w (m/s)",
    title="Artifact B: unsteady flight near inflection\n(screen off, to show what w_min removes)")
for (i, L) in enumerate(labels)
    scatterlines!(ax2, res[L].px, res[L].pw; color=cols[i], markersize=11, label=L)
end
hlines!(ax2, [0.0]; color=:gray, linestyle=:dash)
vlines!(ax2, [6.0]; color=:black, linestyle=:dot)
text!(ax2, 6.6, minimum(filter(isfinite, vcat([res[L].pw for L in labels]...)));
    text="← masked by w_min", fontsize=11, align=(:left, :bottom))
axislegend(ax2; position=:rb)

save(joinpath(OUT, "w_quality_diagnostics.png"), fig)
@info "wrote w_quality_diagnostics.png"
