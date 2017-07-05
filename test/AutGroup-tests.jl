@testset "Automorphisms" begin
   using Nemo
   G = PermutationGroup(4)

   @testset "AutSymbol" begin
      @test_throws MethodError Groups.AutSymbol("a")
      @test_throws MethodError Groups.AutSymbol("a", 1)
      f = Groups.AutSymbol("a", 1, :(a()), v -> v)
      @test isa(f, Groups.GSymbol)
      @test isa(f, Groups.AutSymbol)
      @test isa(Groups.perm_autsymbol(G([1,2,3,4])), Groups.AutSymbol)
      @test isa(Groups.rmul_autsymbol(1,2), Groups.AutSymbol)
      @test isa(Groups.lmul_autsymbol(3,4), Groups.AutSymbol)
      @test isa(Groups.flip_autsymbol(3), Groups.AutSymbol)
   end

   a,b,c,d = Nemo.gens(FreeGroup(4))
   domain = [a,b,c,d]

   @testset "flip_autsymbol correctness" begin
      @test Groups.flip_autsymbol(1)(domain) == [a^-1, b,c,d]
      @test Groups.flip_autsymbol(2)(domain) == [a, b^-1,c,d]
      @test Groups.flip_autsymbol(3)(domain) == [a, b,c^-1,d]
      @test Groups.flip_autsymbol(4)(domain) == [a, b,c,d^-1]
      @test inv(Groups.flip_autsymbol(1))(domain) == [a^-1, b,c,d]
      @test inv(Groups.flip_autsymbol(2))(domain) == [a, b^-1,c,d]
      @test inv(Groups.flip_autsymbol(3))(domain) == [a, b,c^-1,d]
      @test inv(Groups.flip_autsymbol(4))(domain) == [a, b,c,d^-1]
   end

   @testset "perm_autsymbol correctness" begin
      σ = Groups.perm_autsymbol(G([1,2,3,4]))
      @test σ(domain) == domain
      @test inv(σ)(domain) == domain

      σ = Groups.perm_autsymbol(G([2,3,4,1]))
      @test σ(domain) == [b, c, d, a]
      @test inv(σ)(domain) == [d, a, b, c]

      σ = Groups.perm_autsymbol(G([2,1,4,3]))
      @test σ(domain) == [b, a, d, c]
      @test inv(σ)(domain) == [b, a, d, c]

      σ = Groups.perm_autsymbol(G([2,3,1,4]))
      @test σ(domain) == [b,c,a,d]
      @test inv(σ)(domain) == [c,a,b,d]
   end

   @testset "rmul/lmul_autsymbol correctness" begin
      i,j = 1,2
      r = Groups.rmul_autsymbol(i,j)
      l = Groups.lmul_autsymbol(i,j)
      @test r(domain) == [a*b,b,c,d]
      @test inv(r)(domain) == [a*b^-1,b,c,d]
      @test l(domain) == [b*a,b,c,d]
      @test inv(l)(domain) == [b^-1*a,b,c,d]

      i,j = 3,1
      r = Groups.rmul_autsymbol(i,j)
      l = Groups.lmul_autsymbol(i,j)
      @test r(domain) == [a,b,c*a,d]
      @test inv(r)(domain) == [a,b,c*a^-1,d]
      @test l(domain) == [a,b,a*c,d]
      @test inv(l)(domain) == [a,b,a^-1*c,d]


      i,j = 4,3
      r = Groups.rmul_autsymbol(i,j)
      l = Groups.lmul_autsymbol(i,j)
      @test r(domain) == [a,b,c,d*c]
      @test inv(r)(domain) == [a,b,c,d*c^-1]
      @test l(domain) == [a,b,c,c*d]
      @test inv(l)(domain) == [a,b,c,c^-1*d]


      i,j = 2,4
      r = Groups.rmul_autsymbol(i,j)
      l = Groups.lmul_autsymbol(i,j)
      @test r(domain) == [a,b*d,c,d]
      @test inv(r)(domain) == [a,b*d^-1,c,d]
      @test l(domain) == [a,d*b,c,d]
      @test inv(l)(domain) == [a,d^-1*b,c,d]
   end

   @testset "AutGroup/AutGroupElem constructors" begin
      f = Groups.AutSymbol("a", 1, :(a()), v -> v)
      @test isa(AutGroupElem(f), Groups.GWord)
      @test isa(AutGroupElem(f), AutGroupElem)
      @test isa(AutGroup(FreeGroup(3)), Nemo.Group)
      @test isa(AutGroup(FreeGroup(1)), Groups.FPGroup)
      A = AutGroup(FreeGroup(1))
      @test isa(Nemo.gens(A), Vector{AutGroupElem})
      @test length(Nemo.gens(A)) == 1
      A = AutGroup(FreeGroup(1), special=true)
      @test length(Nemo.gens(A)) == 0
      A = AutGroup(FreeGroup(2))
      @test length(Nemo.gens(A)) == 7
      gens = Nemo.gens(A)

      @test isa(A(Groups.rmul_autsymbol(1,2)), AutGroupElem)
      @test A(Groups.rmul_autsymbol(1,2)) in gens

      @test isa(A(Groups.rmul_autsymbol(2,1)), AutGroupElem)
      @test A(Groups.rmul_autsymbol(2,1)) in gens

      @test isa(A(Groups.lmul_autsymbol(1,2)), AutGroupElem)
      @test A(Groups.lmul_autsymbol(1,2)) in gens

      @test isa(A(Groups.lmul_autsymbol(2,1)), AutGroupElem)
      @test A(Groups.lmul_autsymbol(2,1)) in gens

      @test isa(A(Groups.flip_autsymbol(1)), AutGroupElem)
      @test A(Groups.flip_autsymbol(1)) in gens

      @test isa(A(Groups.flip_autsymbol(2)), AutGroupElem)
      @test A(Groups.flip_autsymbol(2)) in gens

      @test isa(A(Groups.perm_autsymbol(PermutationGroup(2)([2,1]))),
         AutGroupElem)
      @test A(Groups.perm_autsymbol(PermutationGroup(2)([2,1]))) in gens
   end

   A = AutGroup(FreeGroup(4))

   @testset "eltary functions" begin

      f = Groups.perm_autsymbol(G([2,3,4,1]))
      @test (Groups.change_pow(f, 2)).pow == 1
      @test (Groups.change_pow(f, -2)).pow == 1
      @test (inv(f)).pow == 1


      f = Groups.perm_autsymbol(G([2,1,4,3]))
      @test isa(inv(f), Groups.AutSymbol)

      @test_throws DomainError f^-1
      @test_throws MethodError f*f

      @test A(f)^-1 == A(inv(f))
   end

   @testset "reductions/arithmetic" begin
      f = Groups.perm_autsymbol(G([2,3,4,1]))

      f² = Groups.r_multiply(A(f), [f], reduced=false)
      @test Groups.simplify_perms!(f²) == false
      @test f²^2 == A()

      a = A(Groups.rmul_autsymbol(1,2))*Groups.flip_autsymbol(2)
      b = Groups.flip_autsymbol(2)*A(inv(Groups.rmul_autsymbol(1,2)))
      @test a*b == b*a
      @test a^3 * b^3 == A()
   end

   @testset "specific Aut(F4) tests" begin
      N = 4
      G = AutGroup(FreeGroup(N))
      S = G.gens
      @test isa(S, Vector{Groups.AutSymbol})
      S = [G(s) for s in unique(S)]
      @test isa(S, Vector{AutGroupElem})
      @test S == Nemo.gens(G)
      @test length(S) == 51
      S_inv = [S..., [inv(s) for s in S]...]
      @test length(unique(S_inv)) == 75

      G = AutGroup(FreeGroup(N), special=true, outer=true)
      S = Nemo.gens(G)
      S_inv = [G(), S..., [inv(s) for s in S]...]
      S_inv = unique(S_inv)
      B_2 = [i*j for (i,j) in Base.product(S_inv, S_inv)]
      @test length(B_2) == 2401
      @test length(unique(B_2)) == 1777
   end

end
