export GeneticCode, StandardGeneticCode, tuple_length, alphabet, inverse

"""
    $(TYPEDEF)

A mutable struct representing a genetic code as a mapping from tuples to labels.
"""
mutable struct GeneticCode{T,U<:BioSequences.NucleicAcid,V<:NucleicAcidAlphabet}
    "The order of the letter of the alphabet."
    alphabet_order::Vector{U}
    "A vector of labels, representing e.g. the assignment of amino acids to all possible codons."
    label_order::NamedArray{T, 1, Vector{T}, Tuple{NamedArrays.OrderedDict{LongSequence{V}, Int64}}}
end

"""
    $(TYPEDSIGNATURES)

A genetic code that uses tuples of length `tuple_length` with 
the alphabet `alphabet`. `label_order` represents the list of labels for each tuple.

The order of `label_order` must be identical to the order of the tuples in 
[`alltuples`](@ref)`(alphabet, tuple_length)`.

```julia
gc = GeneticCode(2, [DNA_A, DNA_T], [AA_A, AA_L, AA_K, AA_P])
```
"""
function GeneticCode(tuple_length::Int, alphabet::Vector, label_order::Vector)
    tuples_list = alltuples(alphabet, tuple_length)
    GeneticCode(alphabet, NamedArray(label_order, (tuples_list,), ("aminoacid",)))
end


"""
    $(TYPEDSIGNATURES)

An empty genetic code filled with identical labels as specified 
with the parameter `init` (default is Stop-Signal *).

```julia
gc1 = GeneticCode(2, [DNA_A, DNA_T]) # only stop signal
gc2 = GeneticCode(2, [DNA_A, DNA_T]; init=AA_L) # only Leucine
```
"""
function GeneticCode(tuple_length::Int, alphabet::Vector; init=AA_Term)
    n = length(alphabet)^tuple_length
    labels = fill(init, n)
    GeneticCode(tuple_length, alphabet, labels)
end

"""
    $(TYPEDSIGNATURES)

Genetic code based on 64 codons according the the NCBI list. 
See `BioSequences.ncbi_trans_table` for a list of all known genetic codes.
The order of the bases is TCAG (or UCAG) (and not alphabetically).
"""
function GeneticCode(code::BioSequences.GeneticCode)
    ucag = [DNA_T, DNA_C, DNA_A, DNA_G]
    GeneticCode(
        3,
        ucag, # [stripped_alphabet(DNA)...],
        map(
            codon -> BioSequences.translate(
                codon;
                code,
                allow_ambiguous_codons=true,
                alternative_start=false,
            )[1],
            alltuples(ucag, 3) # BioCodes.codons
        )
    )
end

"""
    $(TYPEDSIGNATURES)

The Standard Genetic Code. See [`GeneticCode`](@ref) for details.
"""
StandardGeneticCode() = GeneticCode(BioSequences.ncbi_trans_table[1])

Base.:(==)(gc1::GeneticCode, gc2::GeneticCode) = gc1.label_order == gc2.label_order

Base.getindex(gc::GeneticCode, idx) = gc.label_order[idx]

function Base.setindex!(gc::GeneticCode, value, idx)
    gc.label_order[idx] = value
end

"""
    $(TYPEDSIGNATURES)

The number of tuples in the code.

```jldoctest
using BioCodes
length(StandardGeneticCode())
# output
64
```
"""
Base.length(gc::GeneticCode) = length(gc.label_order)

"""
    $(TYPEDSIGNATURES)

The length of the tuples in the code.

```jldoctest
using BioCodes
tuple_length(StandardGeneticCode()) # i.e. codons
# output
3
```
"""
tuple_length(gc::GeneticCode) = length(names(gc.label_order)[1][1])

import BioSymbols.alphabet
"""
    $(TYPEDSIGNATURES)

The alphabet used in the code.

```jldoctest
using BioCodes, BioSymbols
alphabet(StandardGeneticCode())
# output
4-element Vector{DNA}:
 DNA_T
 DNA_C
 DNA_A
 DNA_G
```
"""
BioSymbols.alphabet(gc::GeneticCode) = gc.alphabet_order

import Base.replace!

function Base.replace!(gc::GeneticCode, p::Pair{T,U}) where {T<:LongSequence,U<:Any}
    tuple, label = p
    gc.label_order[tuple] = label
end

"""
    $(TYPEDEF)

A helper structure that stores a single mapping of a genetic code for printing.
I.e. the assignment of a codon to an amino acid.
The data types, however, are not restricted to biological entities.
"""
struct GeneticCodeCell
    tuple
    label
end

function Base.:(==)(gcc1::GeneticCodeCell, gcc2::GeneticCodeCell)
    gcc1.tuple == gcc2.tuple && gcc1.label == gcc2.label
end

Base.show(io::IO, gcc::GeneticCodeCell) = print(io, "$(gcc.tuple): $(gcc.label) ")

"""
    $(TYPEDSIGNATURES)

Converts a genetic code into a matrix where the cells are of type [`GeneticCodeCell`](@ref).

This is intended for printing purposes. Only tuple lengths 2 and 3 are supported yet.
"""
function Base.Matrix(gc::GeneticCode; order::Union{Vector{Int},Nothing}=nothing)
    l = tuple_length(gc) # tuple length
    if (l < 2 || l > 3)
        throw("Unsupported tuple length: $l")
    end

    # The type of the tuples (RNA or DNA sequence):
    T = typeof(names(gc.label_order)[1][1])

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
            return T([Σ[i1], Σ[i2]])
        elseif l == 3
            i1 = (i - 1) ÷ N + 1
            i2 = j
            i3 = ((i - 1) % N) + 1
            return T([Σ[i1], Σ[i2], Σ[i3]])
        end
    end

    function output(i, j)
        t = build_tuple(i, j, l)
        return GeneticCodeCell(t, gc[t]) # "$t: $(gc[t])"
    end

    [output(i, j) for i in 1:nrows, j in 1:ncols]
end

function Base.show(io::IO, gc::GeneticCode)
    M = Matrix(gc)
    n, m = size(M)
    a = join(alphabet(gc), ",")
    t = tuple_length(gc)
    if n < 70 && m < 10
        a = join(gc.alphabet_order, ",")
        println(io, "GeneticCode over alphabet {$a} with tuple length $t")
        for r in eachrow(M)
            print(io, " ")
            for o in r
                print(io, "$o  ")
            end
            println(io)
        end
    else
        println(io, "GeneticCode over alphabet {$a} with tuple length $t(...)")
    end
end

# Base.show(io::IO, gc::Matrix{GeneticCodeCell}) = print(io, gc)


"""
    $(TYPEDSIGNATURES)

Dictionary which maps a label to a set of corresponding tuples
for a given genetic code `gc`. 

# Examples
The Vertebrate Mitochondrial Code (index 2) is used. This code has four
stop codons.
```jldoctest
using BioCodes, BioSequences
gc = GeneticCode(ncbi_trans_table[2])
a2c = inverse(gc)
sort([a2c[AA_Term]...]) # Stop signal (set is sorted)
# output
4-element Vector{LongSequence{DNAAlphabet{4}}}:
 AGA
 AGG
 TAA
 TAG
```   
"""
function inverse(gc::GeneticCode)
    U = typeof(gc.label_order[1])
    tuples = names(gc.label_order)[1]
    T = typeof(tuples[1])

    tuple_table = Dict{U,Set{T}}()

    # Iterate over all possible tuples:
    for i in eachindex(gc.label_order)
        t = tuples[i]
        label = gc.label_order[i]
        # Add the codon to the corresponding amino acid's list in the dictionary:
        if !haskey(tuple_table, label)
            tuple_table[label] = Set{T}()
        end
        push!(tuple_table[label], t)
    end

    return tuple_table
end
