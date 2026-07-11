# Shared mission registry for the GliderADCP.jl examples.
#
# Each entry locates a sea064 deployment and the few mission-specific facts the
# workflow needs; the processing scripts (`currents.jl`, `realtime_onboard.jl`)
# are otherwise mission-agnostic and treat every mission identically. Latitude is
# derived from the navigation at run time, not stored here.
#
# Fields:
#   label   short mission name, used for output filenames ("M38")
#   dir     deployment folder
#   binary  native .ad2cp path, relative to `dir` (every mission has one)
#   prefix  stream counter for segment disambiguation — folders can hold several
#           missions, so streams are matched as "<prefix>.gli.sub" etc. ("38")
#   netcdf  MIDAS netCDF export relative to `dir` if one exists (for a
#           binary-vs-netCDF reader-parity check), else `nothing`

const GDATA = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData"

const MISSIONS = Dict(
    "m37" => (label="M37",
              dir=joinpath(GDATA, "sea064-20221021-norse-janmayen-complete"),
              binary="ad2cp/sea064_M37.ad2cp", prefix="37", netcdf=nothing),
    "m38" => (label="M38",
              dir=joinpath(GDATA, "sea064-20221102-norse-lofoten-complete"),
              binary="ad2cp/102381_sea064_M38/sea064_M38.ad2cp", prefix="38",
              netcdf="ad2cp/102381_sea064_M38/sea064_M38.ad2cp.00000.nc"),
    "m48" => (label="M48",
              dir=joinpath(GDATA, "sea064-20231112-norse-janmayen-complete"),
              binary="ad2cp/sea064_M48.ad2cp", prefix="48",
              netcdf="ad2cp/sea064_M48.ad2cp.00000.nc"),
    "m59" => (label="M59",
              dir=joinpath(GDATA, "sea064-20240720-nesma-passengers-complete"),
              binary="ad2cp/sea064_M59.ad2cp", prefix="59",
              netcdf="ad2cp/sea064_M59.ad2cp.00000.nc"),
)

# Order in which missions run when no CLI argument selects a subset.
const MISSION_ORDER = ["m37", "m38", "m48", "m59"]

selected_missions() = isempty(ARGS) ? MISSION_ORDER : lowercase.(ARGS)
