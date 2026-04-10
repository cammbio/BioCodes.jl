export alltuples, tuples, circshift, stripped_alphabet

stripped_alphabet(_::Type{DNA}) = (DNA_A, DNA_C, DNA_G, DNA_T)
stripped_alphabet(_::Type{RNA}) = (RNA_A, RNA_C, RNA_G, RNA_U)
stripped_alphabet(_::Type{AminoAcid}) = (AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, 
AA_I, AA_L, AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V)

"""
    $(TYPEDSIGNATURES)

Create all possible tuples of length `l` from the given `alphabet`.
The tuples are returned as a vector of `LongDNA` or `LongRNA` sequences,
depending on the type of the alphabet (`DNA` or `RNA`).
"""
function alltuples(alphabet::Vector{T}, l::Int)::Vector where {T<:BioSymbol}
    S = typeof(alphabet[1])
    v = repeat([alphabet], l)
    ts = reshape(collect(Iterators.product(v...)), 1, :)[1, :] #  make vector
    if S == DNA
        return [LongDNA{4}(t) for t in ts]
    else
        return [LongRNA{4}(t) for t in ts]
    end
end

alltuples(alphabet::NTuple{N,T}, l::Int) where {N, T<:BioSymbol} = alltuples([alphabet...], l)

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