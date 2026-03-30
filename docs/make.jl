# Run with: julia --project=./docs/make.jl
using Documenter, DocStringExtensions, BioCodes
makedocs(format=Documenter.HTML(), modules=[BioCodes], sitename="BioCodes.jl",
    authors="Markus Gumbel and other contributors.")
deploydocs(
    repo="github.com/cammbio/BioCodes.jl.git",
    # devbranch = "master",  # or "master", depending on your default branch
    push_preview=true,
    #deps=nothing,
    #make=nothing
)