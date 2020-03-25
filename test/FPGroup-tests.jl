@testset "FPGroups definitions" begin
    F = FreeGroup(["a", "b", "c"])
    a,b,c = Groups.gens(F)
    R = [a^2, a*b*a, c*b*a]
    @test F/R isa FPGroup
    @test F isa FreeGroup
    G = F/R
    A,B,C = Groups.gens(G)

    @test Groups.reduce!(A^2) == one(G)
    @test Groups.reduce!(A*B*A*A) == A
    @test Groups.reduce!(A*A*B*A) == A

    @test Groups.freepreimage(G) == F
    @test Groups.freepreimage(B^2) == b^2

    @test G/[B^2, C*B*C] isa FPGroup
end
