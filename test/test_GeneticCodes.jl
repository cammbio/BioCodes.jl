using BioCodes
using Test, BioSymbols, BioSequences

@testset "Amino acid to codons" begin
    @testset "Standard DNA" begin
        gc = translateAA2Codons()
        @test gc[AA_F] == Set([dna"TTT", dna"TTC"])
        @test gc[AA_Term] == Set([dna"TAA", dna"TAG", dna"TGA"])
    end

    @testset "Standard RNA" begin
        gc = translateAA2Codons(ncbi_trans_table[1], RNA)
        @test gc[AA_F] == Set([rna"UUU", rna"UUC"])
        @test gc[AA_Term] == Set([rna"UAA", rna"UAG", rna"UGA"])
    end
end

@testset "codon to amino acid" begin
    @testset "Standard DNA" begin
        gc = translateCodon2AA()
        @test gc[dna"ATG"] == AA_M
    end
    @testset "Standard RNA" begin
        gc = translateCodon2AA(ncbi_trans_table[1], RNA)
        @test gc[rna"AUG"] == AA_M
    end
end

@testset "GeneticCode" begin
    @testset "Empty 2x2 code" begin
        gc = GeneticCode(2, [DNA_A, DNA_T])
        @test gc.label_order[dna"AT"] == AA_Term
        @test gc[dna"AA"] == AA_Term
        @test length(gc.label_order) == 4
        @test length(gc) == 4
    end

    @testset "SCG" begin
        gc = StandardGeneticCode()
        @test length(gc) == 64
        @test tuple_length(gc) == 3
        @test length(alphabet(gc)) == 4
        @test gc[dna"ATG"] == AA_M

        M = Matrix(gc)
        M[1, 1] = "AAA: K"
    end
end