using Groups.MatrixGroups

@testset "Matrix Groups" begin
    @testset "SL(n, ℤ)" begin
        SL3Z = SpecialLinearGroup{3}(Int8)

        S = gens(SL3Z)
        union!(S, inv.(S))

        _, sizes = Groups.wlmetric_ball(S; radius = 4)

        @test sizes == [13, 121, 883, 5455]

        E(i, j) = SL3Z([A[MatrixGroups.ElementaryMatrix{3}(i, j, Int8(1))]])

        A = alphabet(SL3Z)
        w = E(1, 2)
        r = E(2, 3)^-3
        s = E(1, 3)^2 * E(3, 2)^-1

        S = [w, r, s]
        S = unique([S; inv.(S)])
        _, sizes = Groups.wlmetric_ball(S; radius = 4)
        @test sizes == [7, 33, 141, 561]

        Logging.with_logger(Logging.NullLogger()) do
            @testset "GroupsCore conformance" begin
                test_Group_interface(SL3Z)
                g = SL3Z(rand(1:length(alphabet(SL3Z)), 10))
                h = SL3Z(rand(1:length(alphabet(SL3Z)), 10))

                test_GroupElement_interface(g, h)
            end
        end

        x = w * inv(SL3Z(word(w)[end:end])) * r

        @test length(word(x)) == length(word(r))
        @test size(x) == (3, 3)
        @test eltype(x) == Int8

        @test contains(sprint(show, SL3Z), "SL{3,Int8}")
        @test contains(
            sprint(show, MIME"text/plain"(), SL3Z),
            "special linear group",
        )
        @test contains(sprint(show, MIME"text/plain"(), x), "∈ SL{3,Int8}")
        @test sprint(print, x) isa String

        @test length(word(x)) == 3
    end

    @testset "Sp(6, ℤ)" begin
        Sp6 = MatrixGroups.SymplecticGroup{6}(Int8)

        Logging.with_logger(Logging.NullLogger()) do
            @testset "GroupsCore conformance" begin
                test_Group_interface(Sp6)
                g = Sp6(rand(1:length(alphabet(Sp6)), 10))
                h = Sp6(rand(1:length(alphabet(Sp6)), 10))

                test_GroupElement_interface(g, h)
            end
        end

        x = gens(Sp6, 1) * gens(Sp6, 2)^2
        x *= inv(gens(Sp6, 2)^2) * gens(Sp6, 3)

        @test length(word(x)) == 2
        @test size(x) == (6, 6)
        @test eltype(x) == Int8

        @test contains(sprint(show, Sp6), "Sp{6,Int8}")
        @test contains(
            sprint(show, MIME"text/plain"(), Sp6),
            "group of 6×6 symplectic matrices",
        )
        @test contains(sprint(show, MIME"text/plain"(), x), "∈ Sp{6,Int8}")
        @test sprint(print, x) isa String

        @test length(word(x)) == 2

        for g in gens(Sp6)
            @test MatrixGroups.issymplectic(MatrixGroups.matrix(g))
        end
    end

    @testset "General matrix group" begin
        Sp6 = MatrixGroups.SymplecticGroup{6}(Int8)
        G = Groups.MatrixGroup{6}(Matrix{Int16}.(gens(Sp6)))

        Logging.with_logger(Logging.NullLogger()) do
            @testset "GroupsCore conformance" begin
                test_Group_interface(G)
                g = G(rand(1:length(alphabet(G)), 10))
                h = G(rand(1:length(alphabet(G)), 10))

                test_GroupElement_interface(g, h)
            end
        end

        x = gens(G, 1) * gens(G, 2)^3
        x *= gens(G, 2)^-3

        @test length(word(x)) == 1
        @test size(x) == (6, 6)
        @test eltype(x) == Int16

        @test contains(sprint(show, G), "H ⩽ GL{6,Int16}")
        @test contains(
            sprint(show, MIME"text/plain"(), G),
            "subgroup of 6×6 invertible matrices",
        )
        @test contains(sprint(show, MIME"text/plain"(), x), "∈ H ⩽ GL{6,Int16}")
        @test sprint(print, x) isa String

        @test length(word(x)) == 1
    end
end
