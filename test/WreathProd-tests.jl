@testset "WreathProducts" begin
   S_3 = PermutationGroup(3)
   S_2 = PermutationGroup(2)
   b = perm"(1,2,3)"
   a = perm"(1,2)"

   @testset "Constructors" begin
      @test Groups.WreathProduct(S_2, S_3) isa AbstractAlgebra.Group
      B3 = Groups.WreathProduct(S_2, S_3)
      @test B3 isa Groups.WreathProduct
      @test B3 isa WreathProduct{3, AbstractAlgebra.Generic.PermGroup{Int}, Int}

      aa = Groups.DirectPowerGroupElem((a^0 ,a, a^2))

      @test Groups.WreathProductElem(aa, b) isa AbstractAlgebra.GroupElem
      x = Groups.WreathProductElem(aa, b)
      @test x isa Groups.WreathProductElem
      @test x isa 
         Groups.WreathProductElem{3, AbstractAlgebra.Generic.perm{Int}, Int}

      @test B3.N == Groups.DirectPowerGroup(S_2, 3)
      @test B3.P == S_3

      @test B3(aa, b) == Groups.WreathProductElem(aa, b)
      w = B3(aa, b)
      @test B3(w) == w
      @test B3(b) == Groups.WreathProductElem(B3.N(), b)
      @test B3(aa) == Groups.WreathProductElem(aa, S_3())

      @test B3((a^0 ,a, a^2), b) isa WreathProductElem

      @test B3((a^0 ,a, a^2), b) == B3(aa, b)
   end

   @testset "Types" begin
      B3 = Groups.WreathProduct(S_2, S_3)

      @test elem_type(B3) == Groups.WreathProductElem{3, perm{Int}, Int}

      @test parent_type(typeof(B3())) == Groups.WreathProduct{3, parent_type(typeof(B3.N.group())), Int}

      @test parent(B3()) == Groups.WreathProduct(S_2,S_3)
      @test parent(B3()) == B3
   end

   @testset "Basic operations on WreathProductElem" begin
      aa = Groups.DirectPowerGroupElem((a^0 ,a, a^2))
      B3 = Groups.WreathProduct(S_2, S_3)
      g = B3(aa, b)

      @test g.p == b
      @test g.n == DirectPowerGroupElem(aa.elts)

      h = deepcopy(g)
      @test h == g
      @test !(g === h)

      g = B3(Groups.DirectPowerGroupElem((a ,a, a^2)), g.p)

      @test g.n[1] == parent(g.n[1])(a)
      @test g != h

      @test hash(g) != hash(h)
   end

   @testset "Group arithmetic" begin
      B4 = Groups.WreathProduct(AdditiveGroup(GF(3)), PermutationGroup(4))

      x = B4((0,1,2,0), perm"(1,2,3)(4)")
      @test inv(x) == B4((1,0,2,0), perm"(1,3,2)(4)")

      y = B4((1,0,1,2), perm"(1,4)(2,3)")
      @test inv(y) == B4((1,2,0,2), perm"(1,4)(2,3)")

      @test x*y == B4((0,2,0,2), perm"(1,3,4)(2)")

      @test y*x == B4((1,2,2,2), perm"(1,4,2)(3)")


      @test inv(x)*y == B4((2,1,2,2), perm"(1,2,4)(3)")

      @test y*inv(x) == B4((1,2,1,0), perm"(1,4,3)(2)")
      
      @test (x*y)^6 == ((x*y)^2)^3

   end

   @testset "Iteration" begin
      B3_a = Groups.WreathProduct(AdditiveGroup(GF(3)), S_3)
      @test order(B3_a) == 3^3*6
      @test collect(B3_a) isa Vector{
      WreathProductElem{3, AddGrpElem{AbstractAlgebra.gfelem{Int}}, Int}}

      B3_m = Groups.WreathProduct(MultiplicativeGroup(GF(3)), S_3)
      @test order(B3_m) == 2^3*6
      @test collect(B3_m) isa Vector{
      WreathProductElem{3, MltGrpElem{AbstractAlgebra.gfelem{Int}}, Int}}
      
      @test length(Set([B3_a, B3_m, B3_a])) == 2

      Wr = WreathProduct(PermutationGroup(2),PermutationGroup(4))

      elts = collect(Wr)
      @test elts isa Vector{Groups.WreathProductElem{4, Generic.perm{Int}, Int}}
      @test order(Wr) == 2^4*factorial(4)      

      @test length(elts) == order(Wr)
      @test all([g*inv(g) == Wr() for g in elts])
      @test all(inv(g*h) == inv(h)*inv(g) for g in elts for h in elts)
   end
   
   @testset "Misc" begin
      B3_a = Groups.WreathProduct(AdditiveGroup(GF(3)), S_3)
      @test string(B3_a) == "Wreath Product of The additive group of Finite field F_3 by Permutation group over 3 elements"
      
      @test string(B3_a(perm"(1,3)")) == "([0,0,0]≀(1,3))"
      
      B3_m = Groups.WreathProduct(MultiplicativeGroup(GF(3)), S_3)
      @test string(B3_m) == "Wreath Product of The multiplicative group of Finite field F_3 by Permutation group over 3 elements"
      
      @test string(B3_m(perm"(1,3)")) == "([1,1,1]≀(1,3))"
      
   end

end
