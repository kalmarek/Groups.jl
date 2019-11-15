@testset "FPGroups definitions" begin
    F = FreeGroup(["a", "b", "c"])
    a,b,c = gens(F)
    R = [a^2, a*b*a, c*b*a]
    @test F/R isa FPGroup
    @test F isa FreeGroup
    G = F/R
    A,B,C = gens(G)

    @test A^2 == one(G)
    @test A*B*A*A == A
    @test A*A*B*A == B*A

    @test G/[B^2, C*B*C] isa FPGroup
end
