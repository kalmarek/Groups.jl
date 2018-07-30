@testset "WreathProducts" begin
   S_3 = PermutationGroup(3)
   R, x = PolynomialRing(QQ, "x")
   F, a = NumberField(x^2 + 1, "a")
   b = S_3([2,3,1])

   @testset "Constructors" begin
      @test isa(Groups.WreathProduct(F, S_3), AbstractAlgebra.Group)
      @test isa(Groups.WreathProduct(F, S_3), Groups.WreathProduct)
      @test isa(Groups.WreathProduct(F, S_3), Groups.WreathProduct{AddGrp{Generic.ResField{Generic.Poly{Rational{BigInt}}}}, Int64})

      aa = Groups.DirectProductGroupElem([a^0 ,a, a^2])

      @test isa(Groups.WreathProductElem(aa, b), AbstractAlgebra.GroupElem)
      @test isa(Groups.WreathProductElem(aa, b), Groups.WreathProductElem)
      @test isa(Groups.WreathProductElem(aa, b), Groups.WreathProductElem{AddGrpElem{Generic.ResF{Generic.Poly{Rational{BigInt}}}}, Int64})

      B3 = Groups.WreathProduct(F, S_3)

      @test B3.N == Groups.DirectProductGroup(F, 3)
      @test B3.P == S_3

      @test B3(aa, b) == Groups.WreathProductElem(aa, b)
      @test B3(b) == Groups.WreathProductElem(B3.N(), b)
      @test B3(aa) == Groups.WreathProductElem(aa, S_3())

      @test B3([a^0 ,a, a^2], perm"(1,2,3)") isa WreathProductElem

      @test B3([a^0 ,a, a^2], perm"(1,2,3)") == B3(aa, b)
   end

   @testset "Types" begin
      B3 = Groups.WreathProduct(F, S_3)

      @test elem_type(B3) == Groups.WreathProductElem{AddGrpElem{elem_type(F)}, Int}

      @test parent_type(typeof(B3())) == Groups.WreathProduct{parent_type(typeof(B3.N.group())), Int}

      @test parent(B3()) == Groups.WreathProduct(F,S_3)
      @test parent(B3()) == B3
   end

   @testset "Basic operations on WreathProductElem" begin
      aa = Groups.DirectProductGroupElem([a^0 ,a, a^2])
      B3 = Groups.WreathProduct(F, S_3)
      g = B3(aa, b)

      @test g.p == b
      @test g.n == DirectProductGroupElem(AddGrpElem.(aa.elts))

      h = deepcopy(g)
      @test h == g
      @test !(g === h)

      g.n[1] = parent(g.n[1])(a)

      @test g.n[1] == parent(g.n[1])(a)
      @test g != h

      @test hash(g) != hash(h)

      g.n[1] = a
      @test g.n[1] == parent(g.n[1])(a)
      @test g != h

      @test hash(g) != hash(h)
   end


   @testset "Group arithmetic" begin
      B4 = Groups.WreathProduct(GF(3), PermutationGroup(4))

      x = B4([0,1,2,0], perm"(1,2,3)(4)")
      @test inv(x) == B4([1,0,2,0], perm"(1,3,2)(4)")

      y = B4([1,0,1,2], perm"(1,4)(2,3)")
      @test inv(y) == B4([1,2,0,2], perm"(1,4)(2,3)")

      @test x*y == B4([0,2,0,2], perm"(1,3,4)(2)")

      @test y*x == B4([1,2,2,2], perm"(1,4,2)(3)")


      @test inv(x)*y == B4([2,1,2,2], perm"(1,2,4)(3)")

      @test y*inv(x) == B4([1,2,1,0], perm"(1,4,3)(2)")

   end

   @testset "Misc" begin
      B3 = Groups.WreathProduct(GF(3), S_3)
      @test order(B3) == 3^3*6

      B3 = Groups.WreathProduct(MultiplicativeGroup(GF(3)), S_3)
      @test order(B3) == 2^3*6

      Wr = WreathProduct(PermutationGroup(2),PermutationGroup(4))

      @test isa([elements(Wr)...], Vector{Groups.WreathProductElem{Generic.perm{Int}, Int}})
      @test order(Wr) == 2^4*factorial(4)

      elts = [elements(Wr)...]

      @test length(elts) == order(Wr)
      @test all([g*inv(g) for g in elts] .== Wr())
      @test all(inv(g*h) == inv(h)*inv(g) for g in elts for h in elts)
   end

end
