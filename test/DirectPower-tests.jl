@testset "DirectPowers" begin

   ×(a,b) = Groups.DirectPower(a,b)

   @testset "Constructors" begin
      G = PermutationGroup(3)

      @test Groups.DirectPowerGroup(G,2) isa AbstractAlgebra.Group
      @test G×G isa AbstractAlgebra.Group
      @test Groups.DirectPowerGroup(G,2) isa Groups.DirectPowerGroup{2, Generic.PermGroup{Int64}}

      @test (G×G)×G == DirectPowerGroup(G, 3)
      @test (G×G)×G == (G×G)×G

      GG = DirectPowerGroup(G,2)
      @test (G×G)() isa GroupElem
      @test (G×G)((G(), G())) isa GroupElem
      @test (G×G)([G(), G()]) isa GroupElem

      @test Groups.DirectPowerGroupElem((G(), G())) == (G×G)()
      @test GG(G(), G()) == (G×G)()

      g = perm"(1,2,3)"

      @test GG(g, g^2) isa GroupElem
      @test GG(g, g^2) isa Groups.DirectPowerGroupElem{2, Generic.perm{Int64}}

      h = GG(g,g^2)

      @test h == GG(h)

      @test GG(g, g^2) isa GroupElem
      @test GG(g, g^2) isa Groups.DirectPowerGroupElem

      @test_throws MethodError GG(g,g,g)
      @test GG(g,g^2) == h

      @test h[1] == g
      @test h[2] == g^2
      h = GG(g, G())
      @test h == GG(g, G())
   end

   @testset "Basic arithmetic" begin
      G = PermutationGroup(3)
      GG = G×G
      i = perm"(1,3)"
      g = perm"(1,2,3)"
      
      h = GG(g,g^2)
      k = GG(g^3, g^2)

      @test h^2 == GG(g^2,g)
      @test h^6 == GG()

      @test h*h == h^2
      @test h*k == GG(g,g)

      @test h*inv(h) == (G×G)()
      
      w = GG(g,i)*GG(i,g)
      @test w == GG(perm"(1,2)(3)", perm"(2,3)")
      @test w == inv(w)
      @test w^2 == w*w == GG()
   end

   @testset "elem/parent_types" begin
      G = PermutationGroup(3)
      g = perm"(1,2,3)"

      @test elem_type(G×G) == DirectPowerGroupElem{2, elem_type(G)}
      @test elem_type(G×G×G) == DirectPowerGroupElem{3, elem_type(G)}
      @test parent_type(typeof((G×G)(g,g^2))) == Groups.DirectPowerGroup{2, typeof(G)}
      @test parent(DirectPowerGroupElem((g,g^2,g^3))) == DirectPowerGroup(G,3)
   end

   @testset "Misc" begin
      G = PermutationGroup(3)
      GG = Groups.DirectPowerGroup(G,3)
      @test order(GG) == 216

      @test isa(collect(GG), Vector{Groups.DirectPowerGroupElem{3, elem_type(G)}})
      elts = vec(collect(GG))

      @test length(elts) == 216
      @test all([g*inv(g) == GG() for g in elts])
      @test all(inv(g*h) == inv(h)*inv(g) for g in elts for h in elts)
   end
end
