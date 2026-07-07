# Layer 5 — netCDF export with provenance.

const _REFERENCES = "Visbeck (2002) doi:10.1175/1520-0426(2002)019<0794:DVPULA>2.0.CO;2; " *
    "Todd et al. (2017) doi:10.1175/JTECH-D-16-0156.1; " *
    "Gradone et al. (2023) doi:10.1029/2022JC019608; " *
    "Queste et al., gliderad2cp (JOSS, doi:10.21105/joss.08342)"

"""
    export_sections(path, sections; attrs=Dict{String,Any}())

Write gridded velocity sections (the output of [`grid_profiles`](@ref)) to netCDF with
CF-style metadata and provenance (package version, creation time, references, plus any
user `attrs`). Overwrites `path`.
"""
function export_sections(path::AbstractString, sec; attrs::Dict{String,Any}=Dict{String,Any}())
    NCDataset(path, "c") do ds
        defDim(ds, "depth", length(sec.z))
        defDim(ds, "time", length(sec.t))
        zv = defVar(ds, "depth", Float64.(sec.z), ("depth",))
        zv.attrib["units"] = "m"
        zv.attrib["positive"] = "down"
        tv = defVar(ds, "time", sec.t, ("time",))
        tv.attrib["long_name"] = "segment midpoint time"
        for (name, A, long) in (("u", sec.U, "eastward ocean velocity"),
                                ("v", sec.V, "northward ocean velocity"))
            v = defVar(ds, name, Float64.(A), ("depth", "time"); fillvalue=NaN)
            v.attrib["units"] = "m s-1"
            v.attrib["long_name"] = long
        end
        nv = defVar(ds, "nobs", Int32.(sec.Nobs), ("depth", "time"))
        nv.attrib["long_name"] = "ADCP samples per bin"
        ds.attrib["title"] = "Glider AD2CP absolute ocean velocity sections"
        ds.attrib["source"] = "GliderADCP.jl v$(pkgversion(GliderADCP))"
        ds.attrib["history"] = "created $(Dates.now(UTC))Z"
        ds.attrib["references"] = _REFERENCES
        for (k, v) in attrs
            ds.attrib[k] = v
        end
    end
    return path
end
