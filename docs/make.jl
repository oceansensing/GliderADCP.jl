using Documenter, GliderADCP

makedocs(
    sitename = "GliderADCP.jl",
    modules = [GliderADCP],
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "QA/QC guide" => "qaqc.md",
        "API" => "api.md",
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs],
)
