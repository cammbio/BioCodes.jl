```@meta
# Information for Documenter
CurrentModule = BioCodes
```

```@contents
Pages = ["index.md"]
```

# Tutorial

This tutorial shows how to use `BioCodes.jl`. Since this package is primarily a core package for other packages, most
of its functions are of little practical use. Many features for the analysis of biological codes 
are already available in the [`BioSequences.jl`](https://github.com/BioJulia/BioSequences.jl) package. This tutorial also describes them.

## Working with tuples or codons

### Pre-defined tuples

First, we need to load the required packages:

```@example rt
using BioCodes, BioSequences
```

We can access pre-defined sets of k-mers like codons, dinucleotides, and tetranucleotides:

```@example rt
BioCodes.codons
```

```@example rt
BioCodes.dinucs
```

```@example rt
BioCodes.tetranucs
```

It also possible, to create other sets, e.g. all di-nucleoties over the
alphabet {A, T}:
```@example rt
alltuples((DNA_A, DNA_T), 2)
```

### Tuples from sequences

A sequence can be split into non-overlapping tuples of a given length using the
[`Base.split`](@ref) function. Here, we generate a random DNA sequence of
length 10 and split it into codons (tuples of length 3):

```@example rt
seq = randseq(DNAAlphabet{4}(), 10)
println(seq)
codons = split(seq; l=3)
println(codons)
```
Note that the last codon is incomplete, because the sequence length is not a multiple of 3.

### Shifting sequences

`circshift` shifts a sequence to the left by one position:

```@example rt
circshift(dna"AGCT")
```

Let's shift some RNA in the other direction by explicitly setting the shift direction and step width:

```@example rt
circshift(rna"AGCU"; k=-1)
```

We can also use [`BioCodes.circshift`](@ref) together with [`Base.split`](@ref) to generate all
cyclic permutations of k-mers in a sequence:

```@example rt
seq = rna"CAGCUUGAG"
join(circshift.(split(seq, l=3)))
```

### Complementary and reversed sequences

Complemented tuples can be generated with `BioSequences.complement` which is applied element-wise to a sequence of codons:

```@example rt
seq = rna"CAGCUUGAG"
complement.(split(seq, l=3))
```

Reversed tuples are also available via `BioSequences.reverse`:

```@example rt
seq = rna"CAGCUUGAG"
reverse.(split(seq, l=3))
```

If applied together, we get complementary reversed codons:

```@example rt
seq = rna"CAGCUUGAG"
complement.(reverse.(split(seq, l=3)))
```
Or, if the sequence of codons is supposed to be reserved, too:

```@example rt
seq = rna"CAGCUUGAG"
reverse(complement.(reverse.(split(seq, l=3))))
```

## Genetic code tables

`BioCodes` supports the work with genetic codes which can be known codes with 64 codons 
or artificial codes with an arbitrary number of bases (the alphabet) and tuple sizes.

### Standard genetic code and its derivations

We begin with the standard genetic code:

```@example rt
sgc = StandardGeneticCode()
```

We can also specify the genetic code as listed at NCBI <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>.
`BioSequences.ncbi_trans_table` provides part of it as a list:

```@example rt
using BioSequences

ncbi_trans_table
```

We show the vertebrate mitochondrial code:
```@example rt
GeneticCode(BioSequences.ncbi_trans_table[2])
```

A genetic codes works like an associated array, so it is possible to access the mapping directly:

```@example rt
println(sgc[dna"AAA"])
println(sgc[dna"ATG"]) # start codon
println(sgc[dna"GCG"])
```

Index-based access is also possible. However, this is dependent of the underlying order of the tuples (codons)
and labels (amino acids) - so be aware.

```@example rt
println(sgc[1]) # TTT
println(sgc[64]) # GGG
```

`BioCodes` provides also inverse mappings for genetic code tables. [`BioCodes.inverse`](@ref) creates a dictionary which maps 
an amino acid (or label in general) to a set of associated codons (or tuples in general). This might be useful if a frequent lookup for the mapping is required. As an example a translation is created for the vertebrate mitochondrial code (index 2) which has four stop codons.

```@example rt
gc = GeneticCode(ncbi_trans_table[2])
a2c = inverse(gc)

println(a2c[AA_V])
println(a2c[AA_Term]) # 4 stop codons
```

### Artificial genetic code tables

The standard genetic code as discussed in the previous section is a special case of general genetic codes.
The alphabet and tuple length and even the mapping can be changed.

We start with a code that only uses Adenin and Thymine as the alphabet and that has a tuple length of 2:

```@example rt
gc = GeneticCode(2, [DNA_A, DNA_T])
```
Next we show how this empty code can be modified. We assign new amino acids:

```@example rt
gc[dna"AA"] = AA_F
gc[dna"TT"] = AA_L
gc
```

# API

```@autodocs
Modules = [BioCodes]
Order   = [:type, :function]
```
```@docs
BioCodes.dinucs
```
```@docs
BioCodes.codons
```
```@docs
BioCodes.tetranucs
```