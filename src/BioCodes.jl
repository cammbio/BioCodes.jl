module BioCodes

using BioSequences, BioSymbols
using NamedArrays
using DocStringExtensions
using Pipe: @pipe

# ------------------------------------------------ EXT --------------------------------------------------
# Functions used in the extension:

"""
    Plots.plot(gc::GeneticCode; title=nothing, kwargs...) -> Plots.Plot

Plot a `GeneticCode` as a colored table in the standard textbook layout:

- **Top axis** – second base of the codon/tuple.
- **Left axis** – first base, one label per group of N rows (N = alphabet size).
- **Right axis** – third base, one label per row (only for tuple length 3).
- **Cells** – amino acid (or generic label) with a distinctive background color.

The layout reproduces the appearance of classical genetic code tables in the
literature (e.g. T/U rows × C/A/G columns × T/U/C/A/G inner rows). Only tuple
lengths 2 and 3 are supported.

# Arguments
- `gc`    – a `BioCodes.GeneticCode` object.
- `title` – optional plot title; a descriptive default is generated if omitted.
- Any remaining keyword arguments are forwarded to `Plots.plot`.

# Example
```julia
using BioCodes, BioCodesPlots, Plots
plot(StandardGeneticCode())
plot(StandardGeneticCode(; alphabet=RNAAlphabet))
plot(GeneticCode(2, [DNA_T, DNA_C, DNA_A, DNA_G]))
```
"""
function plot end

include("Tuples.jl")
include("GeneticCodes.jl")

end