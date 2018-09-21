@testset "DirectProducts" begin

   ×(a,b) = Groups.pow(a,b)

   @testset "Constructors" begin
      G = PermutationGroup(3)
      g = G([2,3,1])

      @test Groups.DirectProductGroup(G,2) isa AbstractAlgebra.Group
      @test G×G isa AbstractAlgebra.Group
      @test Groups.DirectProductGroup(G,2) isa Groups.DirectProductGroup{Generic.PermGroup{Int64}}

      @test (G×G)×G == DirectProductGroup(G, 3)
      @test (G×G)×G == (G×G)×G

      F = GF(13)
      FF = F×F
      @test FF×F == F×FF

      GG = DirectProductGroup(G,2)

      @test Groups.DirectProductGroupElem([G(), G()]) == (G×G)()
      @test GG(G(), G()) == (G×G)()

      @test GG([g, g^2]) isa GroupElem
      @test GG([g, g^2]) isa Groups.DirectProductGroupElem{Generic.perm{Int64}}

      h = GG([g,g^2])

      @test h == GG(h)

      @test GG(g, g^2) isa GroupElem
      @test GG(g, g^2) isa Groups.DirectProductGroupElem

      @test_throws DomainError GG(g,g,g)
      @test GG(g,g^2) == h

      @test h[1] == g
      @test h[2] == g^2
      h[2] = G()
      @test h == GG(g, G())

   end

   @testset "Basic arithmetic" begin
      G = PermutationGroup(3)
      g = G([2,3,1])
      h = (G×G)([g,g^2])

      @test h^2 == (G×G)(g^2,g)
      @test h^6 == (G×G)()

      @test h*h == h^2

      @test h*inv(h) == (G×G)()
   end

   @testset "elem/parent_types" begin
      G = PermutationGroup(3)
      g = G([2,3,1])

      @test elem_type(G×G) == DirectProductGroupElem{elem_type(G)}
      @test parent_type(typeof((G×G)(g,g^2))) == Groups.DirectProductGroup{typeof(G)}
      @test parent((G×G)(g,g^2)) == DirectProductGroup(G,2)

      F = AdditiveGroup(GF(13))

      @test elem_type(F×F) == DirectProductGroupElem{Groups.AddGrpElem{AbstractAlgebra.gfelem{Int}}}
      @test parent_type(typeof((F×F)(1,5))) == Groups.DirectProductGroup{Groups.AddGrp{AbstractAlgebra.GFField{Int}}}
      parent((F×F)(1,5)) == DirectProductGroup(F,2)
   end

   @testset "Additive/Multiplicative groups" begin

      R, x = PolynomialRing(QQ, "x")
      F, a = NumberField(x^3 + x + 1, "a")
      G = PermutationGroup(3)

      GG = Groups.DirectProductGroup(G,2)
      FF = Groups.DirectProductGroup(F,2)

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
      FF = Groups.DirectProductGroup(MltGrp(F),2)

      @test FF([a,1]) isa GroupElem
      @test FF([a,1]) isa DirectProductGroupElem
      @test FF([a,1]) isa DirectProductGroupElem{MltGrpElem{elem_type(F)}}
      @test_throws DomainError FF(1,0)
      @test_throws DomainError FF([0,1])
      @test_throws DomainError FF([1,0])

      @test MltGrp(F) isa AbstractAlgebra.Group
      @test MltGrp(F) isa MultiplicativeGroup
      @test DirectProductGroup(MltGrp(F), 2) isa AbstractAlgebra.Group
      @test DirectProductGroup(MltGrp(F), 2) isa DirectProductGroup{MltGrp{typeof(F)}}

      F, a = NumberField(x^3 + x + 1, "a")
      FF = DirectProductGroup(MltGrp(F), 2)

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
      @test DirectProductGroup(AddGrp(F), 2) isa AbstractAlgebra.Group
      @test DirectProductGroup(AddGrp(F), 2) isa DirectProductGroup{AddGrp{typeof(F)}}

      FF = DirectProductGroup(AdditiveGroup(F), 2)

      @test FF([0,a]) isa AbstractAlgebra.GroupElem
      @test FF(F(0),a) isa DirectProductGroupElem
      @test FF(0,0) isa DirectProductGroupElem{AddGrpElem{elem_type(F)}}

      @test FF(F(1),a+1) == FF([1,a+1])

      @test FF([F(1),a+1])*FF([a,a]) == FF(1+a,2a+1)

      x, y = FF([1,a]), FF([a^2,1])
      @test x*y == FF(a^2+1, a+1)
      @test inv(x) == FF([F(-1),-a])

      @test parent(x) == FF
   end

   @testset "Misc" begin
      F = GF(5)

      FF = DirectProductGroup(AdditiveGroup(F),2)
      @test order(FF) == 25

      elts = vec(collect(elements(FF)))
      @test length(elts) == 25
      @test all([g*inv(g) == FF() for g in elts])
      @test all(inv(g*h) == inv(h)*inv(g) for g in elts for h in elts)

      FF = DirectProductGroup(MultiplicativeGroup(F), 3)
      @test order(FF) == 64

      elts = vec(collect(elements(FF)))
      @test length(elts) == 64
      @test all([g*inv(g) == FF() for g in elts])
      @test all(inv(g*h) == inv(h)*inv(g) for g in elts for h in elts)


      G = PermutationGroup(3)
      GG = Groups.DirectProductGroup(G,2)
      @test order(GG) == 36

      @test isa([elements(GG)...], Vector{Groups.DirectProductGroupElem{elem_type(G)}})
      elts = vec(collect(elements(GG)))

      @test length(elts) == 36
      @test all([g*inv(g) == GG() for g in elts])
      @test all(inv(g*h) == inv(h)*inv(g) for g in elts for h in elts)
   end
end
