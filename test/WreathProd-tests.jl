@testset "WreathProducts" begin
   S_3 = PermutationGroup(3)
   F, a = FiniteField(2,3,"a")
   b = S_3([2,3,1])

   @testset "Constructors" begin
      @test isa(Groups.WreathProduct(F, S_3), Nemo.Group)
      @test isa(Groups.WreathProduct(F, S_3), Groups.WreathProduct)
      @test isa(Groups.WreathProduct(F, S_3), Groups.WreathProduct{Nemo.FqNmodFiniteField})

      aa = Groups.DirectProductGroupElem([a^0 ,a, a^2])

      @test isa(Groups.WreathProductElem(aa, b), Nemo.GroupElem)
      @test isa(Groups.WreathProductElem(aa, b), Groups.WreathProductElem)
      @test isa(Groups.WreathProductElem(aa, b), Groups.WreathProductElem{typeof(a)})

      B3 = Groups.WreathProduct(F, S_3)

      @test B3.N == Groups.DirectProductGroup(F, 3)
      @test B3.P == S_3

      @test B3(aa, b) == Groups.WreathProductElem(aa, b)
      @test B3(b) == Groups.WreathProductElem(B3.N(), b)
      @test B3(aa) == Groups.WreathProductElem(aa, S_3())

      g = B3(aa, b)

      @test g.p == b
      @test g.n == aa
      h = deepcopy(g)

      @test hash(g) == hash(h)

      g.n[1] = a

      @test g.n[1] == a
      @test g != h

      @test hash(g) != hash(h)
   end

   @testset "Types" begin
      B3 = Groups.WreathProduct(F, S_3)

      @test elem_type(B3) == Groups.WreathProductElem{elem_type(F), Int}

      @test parent_type(typeof(B3())) == Groups.WreathProduct{parent_type(typeof(B3.N.group())), Int}

      @test parent(B3()) == Groups.WreathProduct(F,S_3)
      @test parent(B3()) == B3
   end

   @testset "Group arithmetic" begin
      B3 = Groups.WreathProduct(F, S_3)

      x = B3(B3.N([1,0,0]), B3.P([2,3,1]))
      y = B3(B3.N([0,1,1]), B3.P([2,1,3]))

      @test x*y == B3(B3.N([0,0,1]), B3.P([3,2,1]))
      @test y*x == B3(B3.N([0,0,1]), B3.P([1,3,2]))

      @test inv(x) == B3(B3.N([0,0,1]), B3.P([3,1,2]))
      @test inv(y) == B3(B3.N([1,0,1]), B3.P([2,1,3]))

      @test inv(x)*y == B3(B3.N([1,1,1]), B3.P([1,3,2]))
      @test y*inv(x) == B3(B3.N([0,1,0]), B3.P([3,2,1]))

   end

   @testset "Misc" begin
      B3 = Groups.WreathProduct(FiniteField(2,1,"a")[1], S_3)
      @test order(B3) == 48

      Wr = WreathProduct(PermutationGroup(2),S_3)

      @test isa([elements(Wr)...], Vector{Groups.WreathProductElem{Generic.perm{Int}, Int}})

      elts = [elements(Wr)...]

      @test length(elts) == order(Wr)
      @test all([g*inv(g) for g in elts] .== Wr())
      @test all(inv(g*h) == inv(h)*inv(g) for g in elts for h in elts)
   end

end
