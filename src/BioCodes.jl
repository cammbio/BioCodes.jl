module BioCodes

export translateAA2Codons, translateCodon2AA, aminoAcidPerCodon

using BioSequences, BioSymbols
using NamedArrays
using DocStringExtensions
using Pipe: @pipe

include("Tuples.jl")
include("GeneticCodes.jl")

end # module end