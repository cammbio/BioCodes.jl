module BioCodes

export translateAA2Codons, translateCodon2AA, aminoAcidPerCodon

using BioSequences, BioSymbols
using NamedArrays
using DocStringExtensions
using Pipe: @pipe

include("Tuples.jl")
include("GeneticCodes.jl")

"""
    $(TYPEDSIGNATURES)

Dictionary which maps an amino acid to a set of corresponding codons 
for a given genetic `code`. The codons can either be `BioSequences.DNA` or 
`BioSequences.RNA` as specified in parameter `S`.
See `BioSequences.ncbi_trans_table` for a list of all known genetic codes.

# Examples
The Vertebrate Mitochondrial Code (index 2) is used. This code has four 
stop codons.
```jldoctest
using BioCodes, BioSequences
a2c = translateAA2Codons(ncbi_trans_table[2], DNA)
sort([a2c[AA_Term]...]) # Stop signal (set is sorted)
# output
4-element Vector{LongSequence{DNAAlphabet{4}}}:
 AGA
 AGG
 TAA
 TAG
```   
"""
function translateAA2Codons(code::BioSequences.GeneticCode, S::Union{Type{DNA},Type{RNA}})
    T = S == DNA ? LongDNA : LongRNA
    codon_table = Dict{AminoAcid,Set{T}}()

    # Iterate over all possible codons (64 in total):
    for codon in alltuples(stripped_alphabet(S), 3)
        # Translate the codon to an amino acid:
        aa = translate(codon; code)[1]
        # Add the codon to the corresponding amino acid's list in the dictionary:
        if !haskey(codon_table, aa)
            codon_table[aa] = Set{T}()
        end
        push!(codon_table[aa], codon)
    end

    return codon_table
end

translateAA2Codons() = translateAA2Codons(ncbi_trans_table[1], DNA)

"""
    $(TYPEDSIGNATURES)

Dictionary which maps a codon to its encoded amino acid
for a given genetic `code`.
See `BioSequences.ncbi_trans_table` for a list of all known genetic codes.

# Examples
The Vertebrate Mitochondrial Code (index 2) is used. This code has four 
stop codons.
```jldoctest
using BioCodes, BioSequences
c2a = translateCodon2AA(ncbi_trans_table[2], DNA)
c2a[dna"ATG"]
# output
AA_M
```   
"""
function translateCodon2AA(code::BioSequences.GeneticCode, S::Union{Type{DNA},Type{RNA}})
    T = S == DNA ? LongDNA : LongRNA
    codon_table = Dict{T,AminoAcid}()

    # Iterate over all possible codons (64 in total):
    for codon in alltuples(stripped_alphabet(S), 3)
        # Translate the codon to an amino acid:
        aa = BioSequences.translate(codon; code)[1]
        push!(codon_table, codon => aa)
    end
    return codon_table
end

translateCodon2AA() = translateCodon2AA(ncbi_trans_table[1], DNA)

end # module end