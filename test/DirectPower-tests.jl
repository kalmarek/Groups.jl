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

      F = AdditiveGroup(GF(13))

      @test elem_type(F×F) == 
         DirectPowerGroupElem{2, Groups.AddGrpElem{AbstractAlgebra.gfelem{Int}}}
      @test parent_type(typeof((F×F)(1,5))) == 
         Groups.DirectPowerGroup{2, Groups.AddGrp{AbstractAlgebra.GFField{Int}}}
      parent((F×F)(1,5)) == DirectPowerGroup(F,2)
      
      F = MultiplicativeGroup(GF(13))

      @test elem_type(F×F) == 
         DirectPowerGroupElem{2, Groups.MltGrpElem{AbstractAlgebra.gfelem{Int}}}
      @test parent_type(typeof((F×F)(1,5))) == 
         Groups.DirectPowerGroup{2, Groups.MltGrp{AbstractAlgebra.GFField{Int}}}
      parent((F×F)(1,5)) == DirectPowerGroup(F,2)
   end

   @testset "Additive/Multiplicative groups" begin

      R, x = PolynomialRing(QQ, "x")
      F, a = NumberField(x^3 + x + 1, "a")
      G = PermutationGroup(3)

      @testset "MltGrp basic functionality" begin
         Gr = MltGrp(F)
         @test Gr(a) isa MltGrpElem
         g = Gr(a)
         @test deepcopy(g) isa MltGrpElem
         @test inv(g) == Gr(a^-1)
         @test Gr() == Gr(1)
         @test inv(g)*g == Gr()
      end

      @testset "AddGrp basic functionality" begin
         Gr = AddGrp(F)
         @test Gr(a) isa AddGrpElem
         g = Gr(a)
         @test deepcopy(g) isa AddGrpElem
         @test inv(g) == Gr(-a)
         @test Gr() == Gr(0)
         @test inv(g)*g == Gr()
      end
   end

   @testset "Direct Product of Multiplicative Groups" begin

      R, x = PolynomialRing(QQ, "x")
      F, a = NumberField(x^3 + x + 1, "a")
      FF = Groups.DirectPowerGroup(MltGrp(F),2)

      @test FF([a,1]) isa GroupElem
      @test FF([a,1]) isa DirectPowerGroupElem
      @test FF([a,1]) isa DirectPowerGroupElem{2, MltGrpElem{elem_type(F)}}
      @test_throws DomainError FF(1,0)
      @test_throws DomainError FF([0,1])
      @test_throws DomainError FF([1,0])

      @test MltGrp(F) isa AbstractAlgebra.Group
      @test MltGrp(F) isa MultiplicativeGroup
      @test DirectPowerGroup(MltGrp(F), 2) isa AbstractAlgebra.Group
      @test DirectPowerGroup(MltGrp(F), 2) isa DirectPowerGroup{2, MltGrp{typeof(F)}}

      F, a = NumberField(x^3 + x + 1, "a")
      FF = DirectPowerGroup(MltGrp(F), 2)

      @test FF(a,a+1) == FF([a,a+1])
      @test FF([1,a+1])*FF([a,a]) == FF(a,a^2+a)
      x, y = FF([1,a]), FF([a^2,1])
      @test x*y == FF([a^2, a])
      @test inv(x) == FF([1,-a^2-1])

      @test parent(x) == FF
   end

   @testset "Direct Product of Additive Groups" begin

      R, x = PolynomialRing(QQ, "x")
      F, a = NumberField(x^3 + x + 1, "a")

      # Additive Group
      @test AddGrp(F) isa AbstractAlgebra.Group
      @test AddGrp(F) isa AdditiveGroup
      @test DirectPowerGroup(AddGrp(F), 2) isa AbstractAlgebra.Group
      @test DirectPowerGroup(AddGrp(F), 2) isa DirectPowerGroup{2, AddGrp{typeof(F)}}

      FF = DirectPowerGroup(AdditiveGroup(F), 2)

      @test FF([0,a]) isa AbstractAlgebra.GroupElem
      @test FF(F(0),a) isa DirectPowerGroupElem
      @test FF(0,0) isa DirectPowerGroupElem{2, AddGrpElem{elem_type(F)}}

      @test FF(F(1),a+1) == FF([1,a+1])

      @test FF([F(1),a+1])*FF([a,a]) == FF(1+a,2a+1)

      x, y = FF([1,a]), FF([a^2,1])
      @test x*y == FF(a^2+1, a+1)
      @test inv(x) == FF([F(-1),-a])

      @test parent(x) == FF
   end

   @testset "Misc" begin
      F = GF(5)

      FF = DirectPowerGroup(AdditiveGroup(F),2)
      @test order(FF) == 25

      elts = vec(collect(FF))
      @test length(elts) == 25
      @test all([g*inv(g) == FF() for g in elts])
      @test all(inv(g*h) == inv(h)*inv(g) for g in elts for h in elts)

      FF = DirectPowerGroup(MultiplicativeGroup(F), 3)
      @test order(FF) == 64

      elts = vec(collect(FF))
      @test length(elts) == 64
      @test all([g*inv(g) == FF() for g in elts])
      @test all(inv(g*h) == inv(h)*inv(g) for g in elts for h in elts)


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
