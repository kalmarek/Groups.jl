@testset "FPGroups" begin
    A = Alphabet([:a, :A, :b, :B, :c, :C], [2,1,4,3,6,5])

    @test New.FreeGroup(A) isa New.FreeGroup
    @test sprint(show, New.FreeGroup(A)) == "free group on 3 generators"

    F = New.FreeGroup([:a, :b, :c], A)
    @test sprint(show, F) == "free group on 3 generators"

    a,b,c = gens(F)
    @test c*b*a isa New.FPGroupElement

    # quotient of F:
    G = New.FPGroup(F, [a*b=>b*a, a*c=>c*a, b*c=>c*b])

    @test G isa New.FPGroup
    @test rand(G) isa New.FPGroupElement

    f = a*c*b
    @test New.word(f) isa Word{UInt8}

    aG,bG,cG = gens(G)

    @test aG isa New.FPGroupElement
    @test_throws AssertionError aG == a
    @test New.word(aG) == New.word(a)

    g = aG*cG*bG

    @test_throws AssertionError f == g
    @test New.word(f) == New.word(g)
    @test New.word(g) == [1, 5, 3]
    New.normalform!(g)
    @test New.word(g) == [1, 3, 5]

    # quotient of G
    H = New.FPGroup(G, [aG^2=>cG, bG*cG=>aG], maxrules=200)

    h = H(New.word(g))

    @test h isa New.FPGroupElement
    @test_throws AssertionError h == g
    @test_throws AssertionError h*g

    New.normalform!(h)
    @test h == H([5])

    @testset "GroupsCore conformance: H" begin
        test_Group_interface(H)
        test_GroupElement_interface(rand(H, 2)...)
    end

end
