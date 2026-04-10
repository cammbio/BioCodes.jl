using BioCodes
using Test, BioSymbols, BioSequences

@testset "Amino acid to codons" begin
    @testset "Standard DNA" begin
        gc = StandardGeneticCode()
        aa2cd = inverse(gc)
        @test aa2cd[AA_F] == Set([dna"TTT", dna"TTC"])
        @test aa2cd[AA_Term] == Set([dna"TAA", dna"TAG", dna"TGA"])
    end
end

@testset "codon to amino acid" begin
    @testset "Standard DNA" begin
        gc = StandardGeneticCode()
        @test gc[dna"ATG"] == AA_M
        @test gc[dna"TGA"] == AA_Term
    end
end

@testset "GeneticCode" begin
    @testset "Empty 2x2 code" begin
        gc = GeneticCode(2, [DNA_A, DNA_T])
        @test gc.label_order[dna"AT"] == AA_Term
        @test length(gc.label_order) == 4
        @test length(gc) == 4

        @test gc[dna"AA"] == AA_Term
        @test gc[dna"AT"] == AA_Term
        @test gc[dna"TA"] == AA_Term
        @test gc[dna"TT"] == AA_Term

        gc[dna"AT"] = AA_K
        @test gc[dna"AA"] == AA_Term
        @test gc[dna"AT"] == AA_K
        @test gc[dna"TA"] == AA_Term
        @test gc[dna"TT"] == AA_Term
    end

    @testset "SCG" begin
        gc = StandardGeneticCode()
        @test length(gc) == 64
        @test tuple_length(gc) == 3
        @test length(alphabet(gc)) == 4
        @test gc[1] == AA_F # as of TCAG order
        @test gc[dna"ATG"] == AA_M

        M = Matrix(gc)
        @test M[1, 1] == BioCodes.GeneticCodeCell(dna"TTT", AA_F)
        @test M[16, 4] == BioCodes.GeneticCodeCell(dna"GGG", AA_G)
    end

    @testset "Nested labels" begin
        gc = GeneticCode(2, [DNA_A, DNA_T]; init=(AA_Term, 42))
        @test length(gc) == 4
        @test tuple_length(gc) == 2
        @test length(alphabet(gc)) == 2
        @test gc[1] == (AA_Term, 42)
        @test gc[dna"TT"] == (AA_Term, 42)

        M = Matrix(gc)
        @test M[1, 1] == BioCodes.GeneticCodeCell(dna"AA", (AA_Term, 42))
        @test M[2, 2] == BioCodes.GeneticCodeCell(dna"TT", (AA_Term, 42))
    end
end