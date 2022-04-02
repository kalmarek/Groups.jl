@testset "GroupConstructions" begin

    @testset "DirectProduct" begin
        GH =
            let G = PermutationGroups.SymmetricGroup(3),
                H = PermutationGroups.SymmetricGroup(4)

                Groups.Constructions.DirectProduct(G, H)
            end
        test_Group_interface(GH)
        test_GroupElement_interface(rand(GH, 2)...)

        @test collect(GH) isa Array{eltype(GH), 2}
        @test contains(sprint(print, GH), "Direct product")
        @test sprint(print, rand(GH)) isa String
    end

    @testset "DirectPower" begin
        GGG = Groups.Constructions.DirectPower{3}(
            PermutationGroups.SymmetricGroup(3),
        )
        test_Group_interface(GGG)
        test_GroupElement_interface(rand(GGG, 2)...)

        @test collect(GGG) isa Array{eltype(GGG), 3}
        @test contains(sprint(print, GGG), "Direct 3-rd power")
        @test sprint(print, rand(GGG)) isa String
    end
    @testset "WreathProduct" begin
        W =
            let G = PermutationGroups.SymmetricGroup(2),
                P = PermutationGroups.SymmetricGroup(4)

                Groups.Constructions.WreathProduct(G, P)
            end
        test_Group_interface(W)
        test_GroupElement_interface(rand(W, 2)...)

        @test collect(W) isa Array{eltype(W), 2}
        @test contains(sprint(print, W), "Wreath product")
        @test sprint(print, rand(W)) isa String
    end
end
