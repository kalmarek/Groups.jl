@testset "Automorphisms" begin

    @testset "Transvections" begin

        @test New.Transvection(:ϱ, 1, 2) isa New.GSymbol
        @test New.Transvection(:ϱ, 1, 2) isa New.Transvection
        @test New.Transvection(:λ, 1, 2) isa New.GSymbol
        @test New.Transvection(:λ, 1, 2) isa New.Transvection
        t = New.Transvection(:ϱ, 1, 2)
        @test inv(t) isa New.GSymbol
        @test inv(t) isa New.Transvection

        @test t != inv(t)

        s = New.Transvection(:ϱ, 1, 2)
        @test t == s
        @test hash(t) == hash(s)

        s_ = New.Transvection(:ϱ, 1, 3)
        @test s_ != s
        @test hash(s_) != hash(s)

        @test New.gersten_alphabet(3) isa Alphabet
        A = New.gersten_alphabet(3)
        @test length(A) == 12

        @test sprint(show, New.ϱ(1, 2)) == "ϱ₁.₂"
        @test sprint(show, New.λ(3, 2)) == "λ₃.₂"
    end

    A4 = Alphabet(
    [:a,:A,:b,:B,:c,:C,:d,:D],
    [ 2, 1, 4, 3, 6, 5, 8, 7]
    )

    A5 = Alphabet(
    [:a,:A,:b,:B,:c,:C,:d,:D,:e,:E],
    [ 2, 1, 4, 3, 6, 5, 8, 7,10, 9]
    )

    F4 = New.FreeGroup([:a, :b, :c, :d], A4)
    a,b,c,d = gens(F4)
    D = ntuple(i->gens(F4, i), 4)

    @testset "Transvection action correctness" begin
        i,j = 1,2
        r = New.Transvection(:ϱ,i,j)
        l = New.Transvection(:λ,i,j)

        (t::New.Transvection)(v::Tuple) = New.evaluate!(v, t, A4)

        @test      r(deepcopy(D)) == (a*b,   b, c, d)
        @test inv(r)(deepcopy(D)) == (a*b^-1,b, c, d)
        @test      l(deepcopy(D)) == (b*a,   b, c, d)
        @test inv(l)(deepcopy(D)) == (b^-1*a,b, c, d)

        i,j = 3,1
        r = New.Transvection(:ϱ,i,j)
        l = New.Transvection(:λ,i,j)
        @test      r(deepcopy(D)) == (a, b, c*a,   d)
        @test inv(r)(deepcopy(D)) == (a, b, c*a^-1,d)
        @test      l(deepcopy(D)) == (a, b, a*c,   d)
        @test inv(l)(deepcopy(D)) == (a, b, a^-1*c,d)

        i,j = 4,3
        r = New.Transvection(:ϱ,i,j)
        l = New.Transvection(:λ,i,j)
        @test      r(deepcopy(D)) == (a, b, c, d*c)
        @test inv(r)(deepcopy(D)) == (a, b, c, d*c^-1)
        @test      l(deepcopy(D)) == (a, b, c, c*d)
        @test inv(l)(deepcopy(D)) == (a, b, c, c^-1*d)

        i,j = 2,4
        r = New.Transvection(:ϱ,i,j)
        l = New.Transvection(:λ,i,j)
        @test      r(deepcopy(D)) == (a, b*d,   c, d)
        @test inv(r)(deepcopy(D)) == (a, b*d^-1,c, d)
        @test      l(deepcopy(D)) == (a, d*b,   c, d)
        @test inv(l)(deepcopy(D)) == (a, d^-1*b,c, d)
    end

    A = New.SpecialAutomorphismGroup(F4, maxrules=1000)

    @testset "AutomorphismGroup constructors" begin
        @test A isa New.AbstractFPGroup
        @test A isa New.AutomorphismGroup
        @test KnuthBendix.alphabet(A) isa Alphabet
        @test New.relations(A) isa Vector{<:Pair}
        @test sprint(show, A) == "automorphism group of free group on 4 generators"
    end

    @testset "Automorphisms: hash and evaluate" begin
        @test New.domain(gens(A, 1)) == D
        g, h = gens(A, 1), gens(A, 8)

        @test New.evaluate(g*h) == New.evaluate(h*g)
        @test (g*h).savedhash == zero(UInt)

        @test sprint(show, typeof(g)) == "Automorphism{Groups.New.FreeGroup{Symbol},…}"

        a = g*h
        b = h*g
        @test hash(a) != zero(UInt)
        @test hash(a) == hash(b)
        @test a.savedhash == b.savedhash

        @test length(unique([a,b])) == 1
        @test length(unique([g*h, h*g])) == 1

        # Not so simple arithmetic: applying starting on the left:
        # ϱ₁₂*ϱ₂₁⁻¹*λ₁₂*ε₂ == σ₂₁₃₄

        g = gens(A, 1)
        x1, x2, x3, x4 = New.domain(g)
        @test New.evaluate(g) == (x1*x2, x2, x3, x4)

        g = g*inv(gens(A, 4)) # ϱ₂₁
        @test New.evaluate(g) == (x1*x2, x1^-1, x3, x4)

        g = g*gens(A, 13)
        @test New.evaluate(g) == (x2, x1^-1, x3, x4)
    end

    @testset "Automorphisms: SAut(F₄)" begin
        N = 4
        G = New.SpecialAutomorphismGroup(New.FreeGroup(N))

        S = gens(G)
        @test S isa Vector{<:New.FPGroupElement{<:New.AutomorphismGroup{<:New.FreeGroup}}}

        @test length(S) == 2*N*(N-1)
        @test length(unique(S)) == length(S)

        S_sym = [S; inv.(S)]
        @test length(S_sym) == length(unique(S_sym))

        pushfirst!(S_sym, one(G))

        B_2 = [i*j for (i,j) in Base.product(S_sym, S_sym)]
        @test length(B_2) == 2401
        @test length(unique(B_2)) == 1777

        @test all(g->isone(inv(g)*g), B_2)
        @test all(g->isone(g*inv(g)), B_2)
    end

    @testset "GroupsCore conformance" begin
        test_Group_interface(A)
        g = A(rand(1:length(KnuthBendix.alphabet(A)), 10))
        h = A(rand(1:length(KnuthBendix.alphabet(A)), 10))

        test_GroupElement_interface(g, h)
    end

end

# using Random
# using GroupsCore
#
# A = New.SpecialAutomorphismGroup(New.FreeGroup(4), maxrules=2000, ordering=KnuthBendix.RecursivePathOrder)
#
# # for seed in 1:1000
# let seed = 68
#     N = 14
#     Random.seed!(seed)
#     g = A(rand(1:length(KnuthBendix.alphabet(A)), N))
#     h = A(rand(1:length(KnuthBendix.alphabet(A)), N))
#     @info "seed=$seed" g h
#     @time isone(g*inv(g))
#     @time isone(inv(g)*g)
#     @info "" length(New.word(New.normalform!(g*inv(g)))) length(New.word(New.normalform!(inv(g)*g)))
#     a = commutator(g, h, g)
#     b = conj(inv(g), h) * conj(conj(g, h), g)
#
#     @info length(New.word(a))
#     @info length(New.word(b))
#
#     w = a*inv(b)
#     @info length(New.word(w))
#     New.normalform!(w)
#     @info length(New.word(w))
#
#
#     #
#     # @time ima = New.evaluate(a)
#     # @time imb = New.evaluate(b)
#     # @info "" a b ima imb
#     # @time a == b
# end
