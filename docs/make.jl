using Documenter, DocumenterCitations
using OceananigansLagrangianFilter, OceananigansLagrangianFilter.Utils

# Set up bibliography 
bib_filepath = joinpath(dirname(@__FILE__), "Lagrangianfilter.bib")
bib = CitationBibliography(bib_filepath, style=:authoryear)

# Set up DocTests (TODO)

# Set up math rendering need to fix this
mathengine = MathJax3(Dict(
    :loader => Dict("load" => ["[tex]/physics"]),
    :tex => Dict(
        "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
        "tags" => "ams",
        "packages" => ["base", "ams", "autoload", "physics"],
    ),
))

format = Documenter.HTML(collapselevel = 1,
                          repolink = "https://loisbaker.github.io/OceananigansLagrangianFilter.jl",
                        # canonical = "https://loisbaker.github.io/OceananigansLagrangianFilter/stable/")
                           mathengine = mathengine)
                        #  size_threshold = 2^20,
                        #  assets = String["assets/citations.css"])

# Literate.jl can do this automatically
# example_pages = [
#     "Geostrophic adjustment online"        => "literated/geostrophic_adjustment_online.md",
#     "Geostrophic adjustment offline"        => "literated/geostrophic_adjustment_offline.md",
#     "Shallow water inertial oscillation offline"       => "literated/shallow_water_inertial_oscillation_offline.md",

# ]

methods_pages = [
    "Eulerian and Lagrangian averaging" => "Methods/Eulerian_Lagrangian_definitions.md",
    "Background: PDEs for Lagrangian filtering" => "Methods/filtering_PDEs.md",
    "Online filtering" => "Methods/online_exponential_filtering.md",
    "Offline filtering" => "Methods/offline_exponential_filtering.md",
    "Filter construction" => "Methods/filter_construction.md",
]
pages = [
    "Home" => "index.md",
    "Quick start" => "quick_start.md",
    "Methods" => methods_pages,
    "Examples" => "examples.md",
    "References" => "references.md",
    "Library" => "library.md", # This page will contain all your functions' docstrings

]

makedocs(; sitename="OceananigansLagrangianFilter.jl", 
           modules = [OceananigansLagrangianFilter],
           remotes=nothing,
           format=format,
           plugins = [bib],
           authors="Lois Baker",
           pages=pages,
           warnonly = [:missing_docs], # Recommended to avoid build failure if docstrings are missing
    )
           