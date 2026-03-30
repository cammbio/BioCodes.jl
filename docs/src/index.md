```@meta
# Information for Documenter
CurrentModule = BioCodes
```

```@contents
Pages = ["index.md"]
```

# Tutorial

This tutorial shows how to use `BioCodes.jl`. This package is mainly a core package for other packages, so most
of the functions are of little use.

## Pre-defined tuples

First, we need to load the required packages:

```@example rt
using BioCodes, BioSequences
```

We can access predefined sets of k-mers like codons, dinucleotides, and tetranucleotides:

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

## Tuples from sequences

A sequence can be split into non-overlapping tuples of a given length using the
[`Base.split`](@ref) function. Here, we generate a random DNA sequence of
length 10 and split it into codons (tuples of length 3):

```@example rt
seq = randseq(DNAAlphabet{4}(), 10)
println(seq)
codons = split(seq; l=3)
println(codons)
```

## Shifting sequences

```@example rt
circshift(dna"AGCT")
```

Let's shift some RNA in the other direction:

```@example rt
circshift(rna"AGCU"; k=-1)
```

We can also use [`BioCodes.circshift`](@ref) together with [`Base.split`](@ref) to generate all
cyclic permutations of k-mers in a sequence:

```@example rt
seq = rna"CAGCUUGAG"
join(circshift.(split(seq, l=3)))
```

## Genetic code tables

`BioCodes` provides mappings for genetic code tables. [`BioCodes.translateAA2Codons`](@ref) creates a dictionary which maps 
an amino acid to a set of associated codons. This might be useful if a frequent lookup for the mapping
is required.

```@example rt
a2c = translateAA2Codons()
println(a2c)

println(a2c[AA_G])
println(a2c[AA_Term])
```

We can also specify the genetic code as listed at NCBI <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>.
`BioSequences.ncbi_trans_table` provides part of it as a list:

```@example rt
using BioSequences

ncbi_trans_table
```
Now a translation is created for the vertebrate mitochondrial code (index 2) which has
four stop codons. `translateAA2Codon` expects two parameters: 1) the genetic code table
and 2) which nucleoties should be used (i.e. `BioSequences.DNA` or `BioSequences.RNA`).

```@example rt
a2c = translateAA2Codons(ncbi_trans_table[2], DNA)

println(a2c[AA_V])
println(a2c[AA_Term]) # 4 stop codons
```

We could also create codons based on RNA.

```@example rt
a2cRNA = translateAA2Codons(ncbi_trans_table[2], RNA)

println(a2cRNA[AA_V])
```

It is also possible to obtain the amino acid which is encoded by a codon. [`BioCodes.translateCodon2AA`](@ref)
creates a dictionary which maps a codon to its amino acid. This is similar to `BioSequences.translate`
with the difference that the amino acid is _not_ a sequence but the amino acid itself (of type `AminoAcid`).

```@example rt
c2a = translateCodon2AA() # Standard genetic code / DNA

println(c2a[dna"ATG"])
```

# API

```@autodocs
Modules = [BioCodes]
Order   = [:type, :function]
```
