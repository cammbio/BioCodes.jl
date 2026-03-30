module BioCodes

export alltuples, tuples, circshift, translateAA2Codons, translateCodon2AA, aminoAcidPerCodon
export stripped_alphabet

using BioSequences, BioSymbols
using DocStringExtensions
using Pipe: @pipe

stripped_alphabet(_::Type{DNA}) = (DNA_A, DNA_T, DNA_C, DNA_G)
stripped_alphabet(_::Type{RNA}) = (RNA_A, RNA_U, RNA_C, RNA_G)
stripped_alphabet(_::Type{AminoAcid}) = (AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, 
AA_I, AA_L, AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V)

"""
    $(TYPEDSIGNATURES)

Create all possible tuples of length `l` from the given `alphabet`.
The tuples are returned as a vector of `LongDNA` or `LongRNA` sequences,
depending on the type of the alphabet (`DNA` or `RNA`).
"""
function alltuples(alphabet::NTuple{N,Union{DNA, RNA}}, l::Int)::Vector{Union{LongDNA,LongRNA}} where {N}
    S = typeof(alphabet[1])
    v = repeat([alphabet], l)
    ts = reshape(collect(Iterators.product(v...)), 1, :)[1, :] #  make vector
    if S == DNA
        return [LongDNA{4}(t) for t in ts]
    else
        return [LongRNA{4}(t) for t in ts]
    end
end

const dinucs = alltuples(stripped_alphabet(DNA), 2)
const codons = alltuples(stripped_alphabet(DNA), 3)
const tetranucs = alltuples(stripped_alphabet(DNA), 4)

import Base.split

"""
    $(TYPEDSIGNATURES)

Create tuples of length `l` by splitting a sequence `seq`.
`seq` must implement the `collect` method.
"""
function Base.split(seq::T; l=3) where {T<:BioSequence}
    S = typeof(seq)
    s = collect(seq) # make vector
    n = length(s) # length of sequence
    e = (n ÷ l) * l # end of sequence as a multiple of l
    @pipe s[1:e] |> reshape(_, l, :) |> eachcol |> map(S, _)
end

import Base.circshift

"""
    $(TYPEDSIGNATURES)

Shift the elements in a sequence for one position to the left (AUGC -> UGCA).
"""
function Base.circshift(seq::T; k=1) where {T<:BioSequence}
    l = length(seq)
    h = k % l
    i = h > 0 ? h : l + h
    if i > 0 # Anything to do at all?
        S = typeof(seq)
        a = collect(seq)
        S([a[(i+1):end]..., a[1:i]...])
    else
        seq
    end
end

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