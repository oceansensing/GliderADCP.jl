# Layer 5 — section gridding of long-format profile tables.

"""
    grid_profiles(prof::DataFrame) -> (t, z, U, V, Nobs)

Assemble per-segment profiles (`yo, t_mid, z, u, v, nobs` — the output schema of
[`solve_inverse`](@ref) / [`solve_shear`](@ref)) into depth × time section matrices,
**matched by depth value** (segments whose profiles start at different depths stay
aligned — unlike the row-index assembly in Slocum-AD2CP). Columns are ordered by
`t_mid`; `Nobs` is 0 where a segment has no bin. `fields` selects the value columns
(default `(:u, :v)`; e.g. `(:w, :w)` for [`solve_w`](@ref) tables, then read `sec.U`).
"""
function grid_profiles(prof::DataFrame; fields::Tuple{Symbol,Symbol}=(:u, :v))
    yos = unique(prof.yo)
    tmid = DateTime[first(prof.t_mid[prof.yo .== y]) for y in yos]
    ord = sortperm(tmid)
    yos, tmid = yos[ord], tmid[ord]
    z = sort(unique(prof.z))
    zrow = Dict(zv => k for (k, zv) in enumerate(z))
    U = fill(NaN, length(z), length(yos))
    V = fill(NaN, length(z), length(yos))
    Nobs = zeros(Int, length(z), length(yos))
    col = Dict(y => j for (j, y) in enumerate(yos))
    for r in eachrow(prof)
        k = zrow[r.z]; j = col[r.yo]
        U[k, j] = getproperty(r, fields[1])
        V[k, j] = getproperty(r, fields[2])
        Nobs[k, j] = r.nobs
    end
    return (t=tmid, z=z, U=U, V=V, Nobs=Nobs)
end
