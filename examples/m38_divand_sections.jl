# DIVAnd-mapped velocity sections from GliderADCP.jl inverse profiles.
#
# Reads the per-yo profiles produced by examples/m38_currents.jl and maps them onto a
# continuous time–depth section with DIVAnd (variational analysis; Barth et al. 2014),
# as a publication-quality alternative to the per-segment heatmap.
#
# Run after m38_currents.jl:
#   JULIA_LOAD_PATH="@:@ocean:@stdlib" julia +1.13 --project=. examples/m38_divand_sections.jl

using GliderADCP, CSV, DataFrames, Dates, Statistics
using DIVAnd, CairoMakie

const OUT = joinpath(@__DIR__, "output")
prof = CSV.read(joinpath(OUT, "M38_profiles_inverse.csv"), DataFrame)
prof = prof[prof.nobs .> 10, :]

t0 = minimum(prof.t_mid)
days = Dates.value.(prof.t_mid .- t0) ./ 86400e3
z = prof.z

# analysis grid: 6-hourly × 10 m over the observed span
di = collect(0:0.25:ceil(maximum(days)))
zi = collect(0:10.0:maximum(z) + 10)
xi = [d for d in di, _ in zi]
yi = [zz for _ in di, zz in zi]
mask = trues(size(xi))
pm = fill(1 / step(range(di[1], di[end]; length=length(di))), size(xi))
pn = fill(1 / 10.0, size(xi))

len = (1.5, 60.0)            # correlation lengths: 1.5 days, 60 m
eps2 = 0.2
sections = Dict{Symbol,Matrix{Float64}}()
for comp in (:u, :v)
    f = Float64.(prof[!, comp])
    fm = mean(f)
    fa, _ = DIVAndrun(mask, (pm, pn), (xi, yi), (days, z), f .- fm, len, eps2)
    sections[comp] = fa .+ fm
end

fig = Figure(size=(1500, 700))
for (row, comp, ttl) in ((1, :u, "U (east) — DIVAnd-mapped inverse"),
                         (2, :v, "V (north) — DIVAnd-mapped inverse"))
    ax = Axis(fig[row, 1]; ylabel="depth (m)", yreversed=true, title=ttl,
        xlabel=row == 2 ? "days since $(Date(t0))" : "")
    hm = heatmap!(ax, di, zi, sections[comp]; colormap=:balance, colorrange=(-0.5, 0.5))
    scatter!(ax, days[1:40:end], z[1:40:end]; markersize=1.5, color=(:black, 0.15))
    row == 1 && Colorbar(fig[1:2, 2], hm; label="velocity (m/s)")
end
save(joinpath(OUT, "M38_UV_divand.png"), fig)
println("saved ", joinpath(OUT, "M38_UV_divand.png"))
