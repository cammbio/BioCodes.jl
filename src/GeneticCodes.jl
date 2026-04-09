export GeneticCode, StandardGeneticCode, tuple_length, alphabet

"""
    $(TYPEDEF)

A mutable struct representing a genetic code as a mapping from RNA codons to amino acids.

# Notes
- Throws an error if `aa_order` does not have length 64.
- The mapping order is determined by the order of codons in [`rna_codon_list`](@ref).
- The struct is mutable, allowing modification of the codon-to-amino acid mapping after creation.

```julia
GeneticCode(aa_order::Vector{AminoAcid})
```
# See also
[`rna_codon_list`](@ref)
"""
mutable struct GeneticCode{U<:BioSequences.BioSymbol}
    #"Tuple length"
    #tuple_length::Int
    "The order of the letter of the alphabet."
    alphabet_order::Vector{U}
    "A vector of labels, representing e.g. the assignment of amino acids to all possible codons."
    label_order::NamedArray #{T,Any,Any} does not work
end

function GeneticCode(tuple_length::Int, alphabet::Vector, label_order::Vector)
    tuples_list = alltuples(alphabet, tuple_length)
    GeneticCode(alphabet, NamedArray(label_order, (tuples_list,), ("aminoacid",)))
end

"Empty genetic code completely filled with stop-signals (*)."
function GeneticCode(tuple_length::Int, alphabet::Vector; init=AA_Term)
    n = length(alphabet)^tuple_length
    labels = fill(init, n)
    GeneticCode(tuple_length, alphabet, labels)
end

Base.:(==)(gc1::GeneticCode, gc2::GeneticCode) = gc1.label_order == gc2.label_order

Base.getindex(gc::GeneticCode, idx) = gc.label_order[idx]

Base.length(gc::GeneticCode) = length(gc.label_order)

tuple_length(gc::GeneticCode) = length(names(gc.label_order)[1][1])

import BioSymbols.alphabet
BioSymbols.alphabet(gc::GeneticCode) = gc.alphabet_order

import Base.replace!

function replace!(gc::GeneticCode, p::Pair{LongSequence,Any})
    tuple, label = p
    gc.label_order[tuple] = label
end

"""
    $(TYPEDSIGNATURES)

Genetic code according the the NCBI list. 
See `BioSequences.ncbi_trans_table` for a list of all known genetic codes.
"""
function GeneticCode(code::BioSequences.GeneticCode)
    GeneticCode(
        3,
        [stripped_alphabet(DNA)...],
        map(
            codon -> BioSequences.translate(
                codon;
                code,
                allow_ambiguous_codons=true,
                alternative_start=false,
            )[1],
            BioCodes.codons
        )
    )
end

StandardGeneticCode() = GeneticCode(BioSequences.ncbi_trans_table[1])

function Base.Matrix(gc::GeneticCode; order::Union{Vector{Int},Nothing}=nothing)
    l = tuple_length(gc) # tuple length
    if (l > 3)
        throw("Unsupported tuple length: $l")
    end

    n = length(gc) # number of tuples
    N = length(alphabet(gc))
    order = isnothing(order) ? (1:N) : order
    Σ = [alphabet(gc)[i] for i ∈ order]

    nrows = n ÷ N
    ncols = N

    function build_tuple(i, j, l)
        if l == 2
            i1 = ((i - 1) % N) + 1
            i2 = j
            return LongDNA{4}([Σ[i1], Σ[i2]])
        elseif l == 3
            i1 = (i - 1) ÷ N + 1
            i2 = j
            i3 = ((i - 1) % N) + 1
            return LongDNA{4}([Σ[i1], Σ[i2], Σ[i3]])
        end
    end

    function output(i, j)
        t = build_tuple(i, j, l)
        return "$t: $(gc[t])"
    end

    [output(i, j) for i in 1:nrows, j in 1:ncols]
end

Base.display(gc::GeneticCode) = display(Matrix(gc))

# ----------------------------------------------
"""
    cmp_custom(a, b, alphabet)

Compare two strings `a` and `b` using a custom alphabet order.

- `alphabet` is a `Vector{Char}` giving the ordering of characters.
- Returns `-1` if `a < b`, `0` if `a == b`, `1` if `a > b` (w.r.t. the custom order).
"""
function cmp_custom(a::AbstractString, b::AbstractString, alphabet::Vector{Char})
    # map each character to its rank
    rank = Dict(c => i for (i, c) in pairs(alphabet))
    default_rank = typemax(Int)  # for chars not in `alphabet`

    # compare character by character
    for (ca, cb) in zip(a, b)
        ra = get(rank, ca, default_rank)
        rb = get(rank, cb, default_rank)
        if ra < rb
            return -1
        elseif ra > rb
            return 1
        end
    end

    # all common prefix equal, shorter string is "smaller"
    if length(a) < length(b)
        return -1
    elseif length(a) > length(b)
        return 1
    else
        return 0
    end
end