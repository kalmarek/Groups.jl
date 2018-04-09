@testset "Automorphisms" begin
   using Nemo
   G = PermutationGroup(Int8(4))

   @testset "AutSymbol" begin
      @test_throws MethodError Groups.AutSymbol("a")
      @test_throws MethodError Groups.AutSymbol("a", 1)
      f = Groups.AutSymbol("a", 1, Groups.FlipAut(2))
      @test isa(f, Groups.GSymbol)
      @test isa(f, Groups.AutSymbol)
      @test isa(Groups.perm_autsymbol(Int8.([1,2,3,4])), Groups.AutSymbol)
      @test isa(Groups.rmul_autsymbol(1,2), Groups.AutSymbol)
      @test isa(Groups.lmul_autsymbol(3,4), Groups.AutSymbol)
      @test isa(Groups.flip_autsymbol(3), Groups.AutSymbol)
   end

   a,b,c,d = Nemo.gens(FreeGroup(4))
   D = NTuple{4,FreeGroupElem}([a,b,c,d])

   @testset "flip_autsymbol correctness" begin
      @test Groups.flip_autsymbol(1)(deepcopy(D)) == (a^-1, b,c,d)
      @test Groups.flip_autsymbol(2)(deepcopy(D)) == (a, b^-1,c,d)
      @test Groups.flip_autsymbol(3)(deepcopy(D)) == (a, b,c^-1,d)
      @test Groups.flip_autsymbol(4)(deepcopy(D)) == (a, b,c,d^-1)
      @test inv(Groups.flip_autsymbol(1))(deepcopy(D)) == (a^-1, b,c,d)
      @test inv(Groups.flip_autsymbol(2))(deepcopy(D)) == (a, b^-1,c,d)
      @test inv(Groups.flip_autsymbol(3))(deepcopy(D)) == (a, b,c^-1,d)
      @test inv(Groups.flip_autsymbol(4))(deepcopy(D)) == (a, b,c,d^-1)
   end

   @testset "perm_autsymbol correctness" begin
      σ = Groups.perm_autsymbol([1,2,3,4])
      @test      σ(deepcopy(D)) == deepcopy(D)
      @test inv(σ)(deepcopy(D)) == deepcopy(D)

      σ = Groups.perm_autsymbol([2,3,4,1])
      @test      σ(deepcopy(D)) == (b, c, d, a)
      @test inv(σ)(deepcopy(D)) == (d, a, b, c)

      σ = Groups.perm_autsymbol([2,1,4,3])
      @test      σ(deepcopy(D)) == (b, a, d, c)
      @test inv(σ)(deepcopy(D)) == (b, a, d, c)

      σ = Groups.perm_autsymbol([2,3,1,4])
      @test      σ(deepcopy(D)) == (b, c, a, d)
      @test inv(σ)(deepcopy(D)) == (c, a, b, d)
   end

   @testset "rmul/lmul_autsymbol correctness" begin
      i,j = 1,2
      r = Groups.rmul_autsymbol(i,j)
      l = Groups.lmul_autsymbol(i,j)
      @test      r(deepcopy(D)) == (a*b,   b, c, d)
      @test inv(r)(deepcopy(D)) == (a*b^-1,b, c, d)
      @test      l(deepcopy(D)) == (b*a,   b, c, d)
      @test inv(l)(deepcopy(D)) == (b^-1*a,b, c, d)

      i,j = 3,1
      r = Groups.rmul_autsymbol(i,j)
      l = Groups.lmul_autsymbol(i,j)
      @test      r(deepcopy(D)) == (a, b, c*a,   d)
      @test inv(r)(deepcopy(D)) == (a, b, c*a^-1,d)
      @test      l(deepcopy(D)) == (a, b, a*c,   d)
      @test inv(l)(deepcopy(D)) == (a, b, a^-1*c,d)

      i,j = 4,3
      r = Groups.rmul_autsymbol(i,j)
      l = Groups.lmul_autsymbol(i,j)
      @test      r(deepcopy(D)) == (a, b, c, d*c)
      @test inv(r)(deepcopy(D)) == (a, b, c, d*c^-1)
      @test      l(deepcopy(D)) == (a, b, c, c*d)
      @test inv(l)(deepcopy(D)) == (a, b, c, c^-1*d)


      i,j = 2,4
      r = Groups.rmul_autsymbol(i,j)
      l = Groups.lmul_autsymbol(i,j)
      @test      r(deepcopy(D)) == (a, b*d,   c, d)
      @test inv(r)(deepcopy(D)) == (a, b*d^-1,c, d)
      @test      l(deepcopy(D)) == (a, d*b,   c, d)
      @test inv(l)(deepcopy(D)) == (a, d^-1*b,c, d)
   end

   @testset "AutGroup/Automorphism constructors" begin
      f = Groups.AutSymbol("a", 1, Groups.FlipAut(1))
      @test isa(Automorphism{3}(f), Groups.GWord)
      @test isa(Automorphism{3}(f), Automorphism)
      @test isa(AutGroup(FreeGroup(3)), Nemo.Group)
      @test isa(AutGroup(FreeGroup(1)), Groups.AbstractFPGroup)
      A = AutGroup(FreeGroup(1))
      @test isa(Nemo.gens(A), Vector{Automorphism{1}})
      @test length(Nemo.gens(A)) == 1
      A = AutGroup(FreeGroup(1), special=true)
      @test length(Nemo.gens(A)) == 0
      A = AutGroup(FreeGroup(2))
      @test length(Nemo.gens(A)) == 7
      gens = Nemo.gens(A)

      @test isa(A(Groups.rmul_autsymbol(1,2)), Automorphism)
      @test A(Groups.rmul_autsymbol(1,2)) in gens

      @test isa(A(Groups.rmul_autsymbol(2,1)), Automorphism)
      @test A(Groups.rmul_autsymbol(2,1)) in gens

      @test isa(A(Groups.lmul_autsymbol(1,2)), Automorphism)
      @test A(Groups.lmul_autsymbol(1,2)) in gens

      @test isa(A(Groups.lmul_autsymbol(2,1)), Automorphism)
      @test A(Groups.lmul_autsymbol(2,1)) in gens

      @test isa(A(Groups.flip_autsymbol(1)), Automorphism)
      @test A(Groups.flip_autsymbol(1)) in gens

      @test isa(A(Groups.flip_autsymbol(2)), Automorphism)
      @test A(Groups.flip_autsymbol(2)) in gens

      @test isa(A(Groups.perm_autsymbol([2,1])), Automorphism)
      @test A(Groups.perm_autsymbol([2,1])) in gens
   end

   A = AutGroup(FreeGroup(4))

   @testset "eltary functions" begin

      f = Groups.perm_autsymbol([2,3,4,1])
      @test (Groups.change_pow(f, 2)).pow == 1
      @test (Groups.change_pow(f, -2)).pow == 1
      @test (inv(f)).pow == 1

      f = Groups.perm_autsymbol([2,1,4,3])
      @test isa(inv(f), Groups.AutSymbol)

      @test_throws DomainError f^-1
      @test_throws MethodError f*f

      @test A(f)^-1 == A(inv(f))
   end

   @testset "reductions/arithmetic" begin
      f = Groups.perm_autsymbol([2,3,4,1])

      f² = Groups.r_multiply(A(f), [f], reduced=false)
      @test Groups.simplifyperms!(f²) == false
      @test f²^2 == A()

      a = A(Groups.rmul_autsymbol(1,2))*Groups.flip_autsymbol(2)
      b = Groups.flip_autsymbol(2)*A(inv(Groups.rmul_autsymbol(1,2)))
      @test a*b == b*a
      @test a^3 * b^3 == A()
      g,h = Nemo.gens(A)[[1,8]] # (g, h) = (ϱ₁₂, ϱ₃₂)

      @test Groups.domain(A) == NTuple{4, FreeGroupElem}(gens(A.objectGroup))

      @test (g*h)(Groups.domain(A)) == (h*g)(Groups.domain(A))
      @test (g*h).savedhash != (h*g).savedhash
      a = g*h
      b = h*g
      @test hash(a) == hash(b)
      @test a.savedhash == b.savedhash
      @test length(unique([a,b])) == 1
      @test length(unique([g*h, h*g])) == 1

      # Not so simple arithmetic: applying starting on the left:
      # ϱ₁₂*ϱ₂₁⁻¹*λ₁₂*ε₂ == σ₂₁₃₄

      g = A(Groups.rmul_autsymbol(1,2))
      x1, x2, x3, x4 = Groups.domain(A)
      @test g(Groups.domain(A)) == (x1*x2, x2, x3, x4)
      g = g*inv(A(Groups.rmul_autsymbol(2,1)))
      @test g(Groups.domain(A)) == (x1*x2, x1^-1, x3, x4)
      g = g*A(Groups.lmul_autsymbol(1,2))
      @test g(Groups.domain(A)) == (x2, x1^-1, x3, x4)
      g = g*A(Groups.flip_autsymbol(2))
      @test g(Groups.domain(A)) == (x2, x1, x3, x4)

      @test g(Groups.domain(A)) == A(Groups.perm_autsymbol([2,1,3,4]))(Groups.domain(A))

      @test g == A(Groups.perm_autsymbol([2,1,3,4]))

      g_im = g(Groups.domain(A))
      @test length(g_im[1]) == 5
      @test length(g_im[2]) == 3
      @test length(g_im[3]) == 1
      @test length(g_im[4]) == 1
      @test length.(Groups.reduce!.(g_im)) == (1,1,1,1)

   end

   @testset "specific Aut(F4) tests" begin
      N = 4
      G = AutGroup(FreeGroup(N))
      S = G.gens
      @test isa(S, Vector{Groups.AutSymbol})
      S = [G(s) for s in unique(S)]
      @test isa(S, Vector{Automorphism{N}})
      @test S == Nemo.gens(G)
      @test length(S) == 51
      S_inv = [S..., [inv(s) for s in S]...]
      @test length(unique(S_inv)) == 75

      G = AutGroup(FreeGroup(N), special=true)
      S = Nemo.gens(G)
      S_inv = [G(), S..., [inv(s) for s in S]...]
      S_inv = unique(S_inv)
      B_2 = [i*j for (i,j) in Base.product(S_inv, S_inv)]
      @test length(B_2) == 2401
      @test length(unique(B_2)) == 1777
   end

end
