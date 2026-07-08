# Task 5: real-time vs delayed-mode products, sea064 M38 (NorSE Lofoten Basin).
#
# Question: how much ocean-velocity accuracy is lost by processing the real-time
# $PNOR ASCII telemetry stream (0.01 m/s velocity quantization, 0.1° attitude, no
# accelerometer, no bottom-track records) instead of the full-resolution .ad2cp
# binary recovered after the mission?  Both sides run the IDENTICAL pipeline —
# same nav (gli.sub), CTD sound speed (pld1.sub), QC, declination, shear-bias
# calibration, DAC — so the AD2CP data source is the only difference.
#
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/m38_realtime_vs_delayed.jl

using GliderADCP
using DataFrames, Dates, Statistics, NaNStatistics
using Printf

const MISSION = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete"
const OUT = joinpath(@__DIR__, "output")
mkpath(OUT)
const LAT0 = 69.5

@info "1/4 load both AD2CP sources + shared nav/CTD"
adcp_d = read_ad2cp(joinpath(MISSION, "ad2cp/102381_sea064_M38/sea064_M38.ad2cp"))
adcp_r = load_pnor(joinpath(MISSION, "delayed/pld1/logs"))
cov_d, cov_r = coverage(adcp_d), coverage(adcp_r)
@info "    delayed (binary):   $(cov_d.n) ens, $(cov_d.t_start) → $(cov_d.t_end)"
@info "    real-time (stream): $(cov_r.n) ens, $(cov_r.t_start) → $(cov_r.t_end)"
# M38: the payload stopped writing the stream on 2022-11-27; the binary adds only
# 750 sparse burst ensembles over the following three months.
nav = load_seaexplorer_nav(joinpath(MISSION, "delayed/nav/logs"))
dac = compute_dac(nav)
pld = load_seaexplorer_pld(joinpath(MISSION, "delayed/pld1/logs"); stream="pld1.sub")
ok = findall(i -> !ismissing(pld.LEGATO_SALINITY[i]) && !ismissing(pld.LEGATO_TEMPERATURE[i]) &&
                  !ismissing(pld.LEGATO_PRESSURE[i]), 1:nrow(pld))
ctd_t = datetime2unix.(pld.time[ok])
c_ctd = soundspeed_from_ctd.(Float64.(pld.LEGATO_SALINITY[ok]),
    Float64.(pld.LEGATO_TEMPERATURE[ok]), Float64.(pld.LEGATO_PRESSURE[ok]), 5.0, LAT0)

@info "2/4 identical pipeline on both sources"
prods = Dict{String,DataFrame}()
for (lab, a, look) in (("delayed", adcp_d, :auto), ("realtime", adcp_r, :down))
    apply_soundspeed!(a, soundspeed_correction(a, ctd_t, c_ctd))
    qc!(a)
    # the stream carries no accelerometer, so look direction must be given explicitly
    p = process_pings(a; lat=LAT0, look=look, declination=magnetic_declination(nav, a.t))
    calibrate_shear_bias!(p)
    prods["$(lab)_inv"] = solve_inverse(p, dac)
    prods["$(lab)_shr"] = solve_shear(p, dac)
    prods["$(lab)_w"] = solve_w(p, dac)
end

@info "3/4 agreement on common (yo, z) bins"
function agreement(a, b, col; nmin=10)
    j = innerjoin(a, b; on=[:yo, :z], makeunique=true)
    c1, c2 = j[!, col], j[!, Symbol(col, :_1)]
    m = (j.nobs .> nmin) .&& (j.nobs_1 .> nmin) .&& isfinite.(c1) .&& isfinite.(c2)
    d = c1[m] .- c2[m]
    (j=j, m=m, n=count(m), r=cor(c1[m], c2[m]), rms=sqrt(mean(d .^ 2)), bias=mean(d))
end
stats = Dict(k => agreement(prods["delayed_$s"], prods["realtime_$s"], col; nmin)
             for (k, s, col, nmin) in (("inv u", "inv", :u, 10), ("inv v", "inv", :v, 10),
                                       ("shr u", "shr", :u, 4), ("shr v", "shr", :v, 4),
                                       ("w", "w", :w, 10)))
for k in ("inv u", "inv v", "shr u", "shr v", "w")
    s = stats[k]
    @printf "    %-6s n=%5d  r=%.4f  rms=%.4f m/s  bias=%+.4f m/s\n" k s.n s.r s.rms s.bias
end

@info "4/4 figure"
try
    @eval using CairoMakie
    fig = Figure(size=(1000, 420))
    su = stats["inv u"]
    ax1 = Axis(fig[1, 1]; xlabel="delayed u (m/s)", ylabel="real-time u (m/s)",
        title=@sprintf("inverse u: r=%.4f, rms=%.1f mm/s", su.r, 1000su.rms), aspect=1)
    scatter!(ax1, su.j.u[su.m], su.j.u_1[su.m]; markersize=2, alpha=0.3)
    ablines!(ax1, 0, 1; color=:black, linestyle=:dash)
    ax2 = Axis(fig[1, 2]; xlabel="rms difference (mm/s)", ylabel="depth (m)",
        yreversed=true, title="real-time − delayed, by depth")
    for (key, color) in (("inv u", :dodgerblue), ("shr u", :darkorange))
        s = stats[key]
        zc, rmsz = Float64[], Float64[]
        for z1 in 0:50:950
            mz = s.m .&& (z1 .<= s.j.z .< z1 + 50)
            count(mz) < 30 && continue
            push!(zc, z1 + 25)
            push!(rmsz, 1000 * sqrt(mean((s.j.u[mz] .- s.j.u_1[mz]) .^ 2)))
        end
        lines!(ax2, rmsz, zc; color, label=key)
    end
    axislegend(ax2; position=:rb)
    save(joinpath(OUT, "M38_realtime_vs_delayed.png"), fig; px_per_unit=2)
    @info "    wrote $(joinpath(OUT, "M38_realtime_vs_delayed.png"))"
catch err
    @warn "figure skipped (CairoMakie not on the load path?)" error = err
end
