using BioCodes
using Test, BioSequences

@testset "Predefined set" begin
    @test length(BioCodes.dinucs) == 16
    @test length(BioCodes.codons) == 64
    @test length(BioCodes.tetranucs) == 256

    S = alltuples((DNA_A, DNA_T), 2)
    @test length(S) == 4
end

@testset "split" begin

    @testset "Empty sequence" begin
        c = split(dna"")
        @test length(c) == 0
    end

    @testset "Codons normal" begin
        c = split(dna"ATACGC")
        @test length(c) == 2
        @test c[1] == dna"ATA"
        @test c[2] == dna"CGC"
    end

    @testset "Codons no multiple of 3" begin
        c = split(dna"ATACG")
        @test length(c) == 1
        @test c[1] == dna"ATA"
    end

    @testset "2-tuples no multiple of 2" begin
        c = split(dna"ATACG"; l=2)
        @test length(c) == 2
        @test c[1] == dna"AT"
        @test c[2] == dna"AC"
    end

    @testset "Amino acids" begin
        c = split(aa"MATMCGYUI"; l=4)
        @test length(c) == 2
        @test c[1] == aa"MATM"
        @test c[2] == aa"CGYU"
    end

    # @testset "Strings" begin
    #     c = split("Markus Gumbel")
    #     @test length(c) == 4
    #     @test c[1] == "Mar"
    #     @test c[2] == "kus"
    #     @test c[3] == " Gu"
    #     @test c[4] == "mbe"
    # end
end

@testset "Tuple diffs" begin
    @testset "Equal" begin
        d = dna"AAA" .!= dna"AAA"
        @test d == [false, false, false]
    end

    @testset "One diff" begin
        d = dna"AAA" .!= dna"AAT"
        @test d == [false, false, true]
    end

    @testset "All diff" begin
        d = dna"GTA" .!= dna"CAT"
        @test d == [true, true, true]
    end
end

@testset "circshift" begin
    @test circshift(dna"A") == dna"A"
    @test circshift(dna"AT") == dna"TA"
    @test circshift(dna"ATG") == dna"TGA"
    @test circshift(dna"ATGC") == dna"TGCA"

    @test circshift(dna"A"; k=1) == dna"A"
    @test circshift(dna"AT"; k=1) == dna"TA"
    @test circshift(dna"ATG"; k=1) == dna"TGA"
    @test circshift(dna"ATGC"; k=1) == dna"TGCA"

    @test circshift(dna"A"; k=2) == dna"A"
    @test circshift(dna"AT"; k=2) == dna"AT"
    @test circshift(dna"ATG"; k=2) == dna"GAT"
    @test circshift(dna"ATGC"; k=2) == dna"GCAT"

    @test circshift(dna"A"; k=-1) == dna"A"
    @test circshift(dna"AT"; k=-1) == dna"TA"
    @test circshift(dna"ATG"; k=-1) == dna"GAT"
    @test circshift(dna"ATGC"; k=-1) == dna"CATG"

    @test circshift(dna"A"; k=-2) == dna"A"
    @test circshift(dna"AT"; k=-2) == dna"AT"
    @test circshift(dna"ATG"; k=-2) == dna"TGA"
    @test circshift(dna"ATGC"; k=-2) == dna"GCAT"

    @test circshift(dna"A"; k=-(1 + 1)) == dna"A"
    @test circshift(dna"AT"; k=-(1 + 2)) == dna"TA"
    @test circshift(dna"ATG"; k=-(1 + 3)) == dna"GAT"
    @test circshift(dna"ATGC"; k=-(1 + 4)) == dna"CATG"

    @test circshift(dna"A"; k=-(2 + 1)) == dna"A"
    @test circshift(dna"AT"; k=-(2 + 2)) == dna"AT"
    @test circshift(dna"ATG"; k=-(2 + 3)) == dna"TGA"
    @test circshift(dna"ATGC"; k=-(2 + 4)) == dna"GCAT"
end

@testset "stripped alphabet" begin
    @testset "DNA" begin
        nts = stripped_alphabet(DNA)
        @test length(nts) == 4
        @test DNA_A ∈ nts
        @test !(DNA_B ∈ nts)
    end

    @testset "RNA" begin
        nts = stripped_alphabet(RNA)
        @test length(nts) == 4
        @test RNA_A ∈ nts
        @test !(RNA_B ∈ nts)
    end

    @testset "Amino acids" begin
        aas = stripped_alphabet(AminoAcid)
        @test length(aas) == 20
        @test AA_R ∈ aas
        @test !(AA_Term ∈ aas)
    end
end