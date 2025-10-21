using Documenter, DocumenterCitations
using Oceananigans
using OceananigansLagrangianFilter, OceananigansLagrangianFilter.Utils
using Literate
using Printf

# Set up bibliography 
bib_filepath = joinpath(dirname(@__FILE__), "Lagrangianfilter.bib")
bib = CitationBibliography(bib_filepath, style=:authoryear)

#####
##### Generate examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")

# The examples that take longer to run should be first. This ensures that the
# docs built which extra workers is as efficient as possible.
example_scripts = [
    "online_filter_geostrophic_adjustment.jl",
    "offline_filter_geostrophic_adjustment.jl",
#    "offline_filter_shallow_water_IO.jl",
]


@info string("Executing the examples")

# for n in 1:length(example_scripts)
#     example = example_scripts[n]
#     example_filepath = joinpath(EXAMPLES_DIR, example)
    
#     start_time = time_ns()
#     Literate.markdown(example_filepath, OUTPUT_DIR;
#                         flavor = Literate.DocumenterFlavor(), execute = true)
#     elapsed = 1e-9 * (time_ns() - start_time)
#     @info @sprintf("%s example took %s to build.", example, prettytime(elapsed))
    
# end


# Set up math rendering 
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
example_pages = [
    "Geostrophic adjustment online"        => "literated/online_filter_geostrophic_adjustment.md",
    "Geostrophic adjustment offline"        => "literated/offline_filter_geostrophic_adjustment.md",
  #  "Shallow water inertial oscillation offline"       => "literated/offline_filter_shallow_water_IO.md",
]

theory_pages = [
    "Eulerian and Lagrangian averaging" => "theory/Eulerian_Lagrangian_definitions.md",
    "Background: PDEs for Lagrangian filtering" => "theory/filtering_PDEs.md",
    "Online Lagrangian filtering equations" => "theory/online_equations.md",
    "Offline Lagrangian filtering equations" => "theory/offline_equations.md"
]
online_pages = [
    "Online filtering implementation" => "online_filtering/online_implementation.md",
    "Choosing online filters" => "online_filtering/choosing_online_filters.md",
]
offline_pages = [
    "Offline filtering implementation" => "offline_filtering/offline_implementation.md",
    "Choosing offline filters" => "offline_filtering/choosing_offline_filters.md",
    "How it works" => "offline_filtering/offline_how_it_works.md",
]
pages = [
    "Home" => "index.md",
    "Installation" => "installation.md",
    "Quick start" => "quick_start.md",
    "Theory" => theory_pages,
    "Online filtering" => online_pages,
    "Offline filtering" => offline_pages,
#    "Examples" => example_pages,
    "References" => "references.md",
    "Library" => "library.md", # This page will contain all your functions' docstrings
]

makedocs(; sitename="OceananigansLagrangianFilter.jl", 
           modules = [OceananigansLagrangianFilter],
           remotes=nothing,
           format=format,
           plugins = [bib],
           authors="Lois Baker",
           doctest = true,
           pages=pages,
           warnonly = [:missing_docs], # Recommended to avoid build failure if docstrings are missing
    )
      
deploydocs(
    repo = "github.com/loisbaker/OceananigansLagrangianFilter.jl.git",
)