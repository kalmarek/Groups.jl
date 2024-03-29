@testset "FPGroups" begin
    A = Alphabet([:a, :A, :b, :B, :c, :C], [2, 1, 4, 3, 6, 5])

    @test FreeGroup(A) isa FreeGroup
    @test sprint(show, FreeGroup(A)) == "free group on 3 generators"

    F = FreeGroup([:a, :b, :c], Groups.KnuthBendix.LenLex(A))
    @test sprint(show, F) == "free group on 3 generators"

    a, b, c = gens(F)
    @test c * b * a isa FPGroupElement

    # quotient of F:
    G = FPGroup(F, [a * b => b * a, a * c => c * a, b * c => c * b])

    @test G isa FPGroup
    @test sprint(show, G) ==
          "⟨ a  b  c | \n\t  a*b => b*a  a*c => c*a  b*c => c*b ⟩"
    @test rand(G) isa FPGroupElement

    f = a * c * b
    @test word(f) isa Word{UInt8}

    aG, bG, cG = gens(G)

    @test aG isa FPGroupElement
    @test_throws AssertionError aG == a
    @test word(aG) == word(a)

    g = aG * cG * bG

    @test_throws AssertionError f == g
    @test word(f) == word(g)
    @test word(g) == [1, 5, 3]
    Groups.normalform!(g)
    @test word(g) == [1, 3, 5]

    let g = aG * cG * bG
        # test that we normalize g before printing
        @test sprint(show, g) == "a*b*c"
    end

    # quotient of G
    H = FPGroup(G, [aG^2 => cG, bG * cG => aG]; max_rules = 200)

    h = H(word(g))

    @test h isa FPGroupElement
    @test_throws AssertionError h == g
    @test_throws MethodError h * g

    H′ = FPGroup(G, [aG^2 => cG, bG * cG => aG]; max_rules = 200)
    @test_throws AssertionError one(H) == one(H′)

    Groups.normalform!(h)
    @test h == H([5])

    @test_logs (
        :warn,
        "using generic isfiniteorder(::AbstractFPGroupElement): the returned `false` might be wrong",
    ) isfiniteorder(h)

    @test_logs (
        :warn,
        "using generic isfinite(::AbstractFPGroup): the returned `false` might be wrong",
    ) isfinite(H)

    Logging.with_logger(Logging.NullLogger()) do
        @testset "GroupsCore conformance: H" begin
            test_Group_interface(H)
            test_GroupElement_interface(rand(H, 2)...)
        end
    end

    @testset "hash/normalform #28" begin
        function cyclic_group(n::Integer)
            A = Alphabet([:a, :A], [2, 1])
            F = FreeGroup(A)
            a, = Groups.gens(F)
            e = one(F)
            Cₙ = FPGroup(F, [a^n => e])

            return Cₙ
        end

        n = 15
        G = cyclic_group(n)
        ball, sizes = Groups.wlmetric_ball(gens(G); radius = n)
        @test first(sizes) == 2
        @test last(sizes) == n

        @test Set(ball) == Set([first(gens(G))^i for i in 0:n-1])
    end
end
