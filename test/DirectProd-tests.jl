@testset "DirectProducts" begin
   using Nemo

   G = PermutationGroup(3)
   g = G([2,3,1])
   F, a = FiniteField(2,3,"a")

   @testset "Constructors" begin
      @test isa(Groups.DirectProductGroup(G,2), Nemo.Group)
      @test isa(G×G, Nemo.Group)
      @test isa(Groups.DirectProductGroup(G,2), Groups.DirectProductGroup{Generic.PermGroup{Int64}})

      GG = Groups.DirectProductGroup(G,2)

      @test GG == Groups.DirectProductGroup(G,2)

      @test Groups.DirectProductGroupElem([G(), G()]) == GG()
      @test GG(G(), G()) == GG()

      @test isa(GG([g, g^2]), GroupElem)
      @test isa(GG([g, g^2]), Groups.DirectProductGroupElem{Generic.perm{Int64}})

      h = GG([g,g^2])

      @test h == GG(h)

      @test isa(GG(g, g^2), GroupElem)
      @test isa(GG(g, g^2), Groups.DirectProductGroupElem)

      @test_throws String GG(g,g,g)
      @test GG(g,g^2) == h

      @test size(h) == (2,)
      @test h[1] == g
      @test h[2] == g^2
      h[2] = G()
      @test h == GG(g, G())
   end

   GG = Groups.DirectProductGroup(G,2)
   FF = Groups.DirectProductGroup(F,2)

   @testset "Types" begin
      @test elem_type(GG) == Groups.DirectProductGroupElem{elem_type(G)}
      @test elem_type(FF) == Groups.DirectProductGroupElem{elem_type(F)}
      @test parent_type(typeof(GG(g,g^2))) == Groups.DirectProductGroup{typeof(G)}
      @test parent_type(typeof(FF(a,a^2))) == Groups.DirectProductGroup{typeof(F)}

      @test isa(FF([0,1]), GroupElem)
      @test isa(FF([0,1]), Groups.DirectProductGroupElem)
      @test isa(FF([0,1]), Groups.DirectProductGroupElem{elem_type(F)})
      @test_throws MethodError FF(1,0)
   end

   @testset "Group arithmetic" begin
      g = G([2,3,1])
      h = GG([g,g^2])

      @test h^2 == GG(g^2,g)
      @test h^6 == GG()

      @test h*h == h^2

      @test h*inv(h) == GG()

      @test FF([0,a])*FF([a,1]) == FF(a,1+a)
      x, y = FF([1,a]), FF([a^2,1])
      @test x*y == FF([a^2+1, a+1])
      @test inv(x) == FF([1,a])
   end

   @testset "Misc" begin
      @test order(GG) == 36
      @test order(FF) == 64

      @test isa([elements(GG)...], Vector{Groups.DirectProductGroupElem{elem_type(G)}})
      elts = [elements(GG)...]

      @test length(elts) == 36
      @test all([g*inv(g) for g in elts] .== GG())
      @test all(inv(g*h) == inv(h)*inv(g) for g in elts for h in elts)
   end
end
