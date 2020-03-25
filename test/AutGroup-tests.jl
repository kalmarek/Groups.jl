@testset "Automorphisms" begin

   G = PermutationGroup(Int8(4))

   @testset "AutSymbol" begin
      @test_throws MethodError Groups.AutSymbol(:a)
      @test_throws MethodError Groups.AutSymbol(:a, 1)
      f = Groups.AutSymbol(:a, 1, Groups.FlipAut(2))
      @test isa(f, Groups.GSymbol)
      @test isa(f, Groups.AutSymbol)
      @test isa(Groups.AutSymbol(perm"(4)"), Groups.AutSymbol)
      @test isa(Groups.AutSymbol([2,3,4,1]), Groups.AutSymbol)
      @test isa(Groups.transvection_R(1,2), Groups.AutSymbol)
      @test isa(Groups.transvection_R(3,4), Groups.AutSymbol)
      @test isa(Groups.flip(3), Groups.AutSymbol)
   end

   a,b,c,d = gens(FreeGroup(4))
   D = NTuple{4,FreeGroupElem}([a,b,c,d])

   @testset "flip correctness" begin
      @test Groups.flip(1)(deepcopy(D)) == (a^-1, b,c,d)
      @test Groups.flip(2)(deepcopy(D)) == (a, b^-1,c,d)
      @test Groups.flip(3)(deepcopy(D)) == (a, b,c^-1,d)
      @test Groups.flip(4)(deepcopy(D)) == (a, b,c,d^-1)
      @test inv(Groups.flip(1))(deepcopy(D)) == (a^-1, b,c,d)
      @test inv(Groups.flip(2))(deepcopy(D)) == (a, b^-1,c,d)
      @test inv(Groups.flip(3))(deepcopy(D)) == (a, b,c^-1,d)
      @test inv(Groups.flip(4))(deepcopy(D)) == (a, b,c,d^-1)
   end

   @testset "perm correctness" begin
      σ = Groups.AutSymbol(perm"(4)")
      @test      σ(deepcopy(D)) == deepcopy(D)
      @test inv(σ)(deepcopy(D)) == deepcopy(D)

      σ = Groups.AutSymbol(perm"(1,2,3,4)")
      @test      σ(deepcopy(D)) == (b, c, d, a)
      @test inv(σ)(deepcopy(D)) == (d, a, b, c)

      σ = Groups.AutSymbol(perm"(1,2)(4,3)")
      @test      σ(deepcopy(D)) == (b, a, d, c)
      @test inv(σ)(deepcopy(D)) == (b, a, d, c)

      σ = Groups.AutSymbol(perm"(1,2,3)(4)")
      @test      σ(deepcopy(D)) == (b, c, a, d)
      @test inv(σ)(deepcopy(D)) == (c, a, b, d)
   end

   @testset "rmul/transvection_R correctness" begin
      i,j = 1,2
      r = Groups.transvection_R(i,j)
      l = Groups.transvection_L(i,j)
      @test      r(deepcopy(D)) == (a*b,   b, c, d)
      @test inv(r)(deepcopy(D)) == (a*b^-1,b, c, d)
      @test      l(deepcopy(D)) == (b*a,   b, c, d)
      @test inv(l)(deepcopy(D)) == (b^-1*a,b, c, d)

      i,j = 3,1
      r = Groups.transvection_R(i,j)
      l = Groups.transvection_L(i,j)
      @test      r(deepcopy(D)) == (a, b, c*a,   d)
      @test inv(r)(deepcopy(D)) == (a, b, c*a^-1,d)
      @test      l(deepcopy(D)) == (a, b, a*c,   d)
      @test inv(l)(deepcopy(D)) == (a, b, a^-1*c,d)

      i,j = 4,3
      r = Groups.transvection_R(i,j)
      l = Groups.transvection_L(i,j)
      @test      r(deepcopy(D)) == (a, b, c, d*c)
      @test inv(r)(deepcopy(D)) == (a, b, c, d*c^-1)
      @test      l(deepcopy(D)) == (a, b, c, c*d)
      @test inv(l)(deepcopy(D)) == (a, b, c, c^-1*d)

      i,j = 2,4
      r = Groups.transvection_R(i,j)
      l = Groups.transvection_L(i,j)
      @test      r(deepcopy(D)) == (a, b*d,   c, d)
      @test inv(r)(deepcopy(D)) == (a, b*d^-1,c, d)
      @test      l(deepcopy(D)) == (a, d*b,   c, d)
      @test inv(l)(deepcopy(D)) == (a, d^-1*b,c, d)
   end

   @testset "AutGroup/Automorphism constructors" begin

      f = Groups.AutSymbol(:a, 1, Groups.FlipAut(1))
      @test isa(Automorphism{3}(f), Groups.GWord)
      @test isa(Automorphism{3}(f), Automorphism)
      @test isa(AutGroup(FreeGroup(3)), AbstractAlgebra.Group)
      @test isa(AutGroup(FreeGroup(1)), Groups.AbstractFPGroup)
      A = AutGroup(FreeGroup(1))
      @test isa(Groups.gens(A), Vector{Automorphism{1}})
      @test length(Groups.gens(A)) == 1
      A = AutGroup(FreeGroup(1), special=true)
      @test length(Groups.gens(A)) == 0
      A = AutGroup(FreeGroup(2))
      @test length(Groups.gens(A)) == 7
      Agens = Groups.gens(A)

      @test isa(A(Groups.transvection_R(1,2)), Automorphism)
      @test A(Groups.transvection_R(1,2)) in Agens

      @test isa(A(Groups.transvection_R(2,1)), Automorphism)
      @test A(Groups.transvection_R(2,1)) in Agens

      @test isa(A(Groups.transvection_R(1,2)), Automorphism)
      @test A(Groups.transvection_R(1,2)) in Agens

      @test isa(A(Groups.transvection_R(2,1)), Automorphism)
      @test A(Groups.transvection_R(2,1)) in Agens

      @test isa(A(Groups.flip(1)), Automorphism)
      @test A(Groups.flip(1)) in Agens

      @test isa(A(Groups.flip(2)), Automorphism)
      @test A(Groups.flip(2)) in Agens

      @test isa(A(Groups.AutSymbol(perm"(1,2)")), Automorphism)
      @test A(Groups.AutSymbol(perm"(1,2)")) in Agens
   end

   A = AutGroup(FreeGroup(4))

   @testset "eltary functions" begin

      f = Groups.AutSymbol(perm"(1,2,3,4)")
      @test (Groups.change_pow(f, 2)).pow == 1
      @test (Groups.change_pow(f, -2)).pow == 1
      @test (inv(f)).pow == 1

      f = Groups.AutSymbol(perm"(1,2)(3,4)")
      @test isa(inv(f), Groups.AutSymbol)

      @test_throws MethodError f*f

      @test A(f)^-1 == A(inv(f))
   end

   @testset "reductions/arithmetic" begin
      f = Groups.AutSymbol(perm"(1,2,3,4)")

      f² = append!(A(f), [f])
      @test Groups.simplifyperms!(Bool, f²) == false
      @test f²^2 == one(A)
      @test !isone(f²)

      a = A(Groups.transvection_L(1,2))*Groups.flip(2)
      b = Groups.flip(2)*A(inv(Groups.transvection_L(1,2)))
      @test a*b == b*a
      @test a^3 * b^3 == one(A)
      g,h = Groups.gens(A)[[1,8]] # (g, h) = (ϱ₁₂, ϱ₃₂)

      @test Groups.domain(A) == NTuple{4, FreeGroupElem}(gens(A.objectGroup))

      @test (g*h)(Groups.domain(A)) == (h*g)(Groups.domain(A))
      @test (g*h).savedhash == zero(UInt)
      @test (h*g).savedhash == zero(UInt)
      a = g*h
      b = h*g
      @test hash(a) != zero(UInt)
      @test hash(b) == hash(a)
      @test a.savedhash == b.savedhash
      @test length(unique([a,b])) == 1
      @test length(unique([g*h, h*g])) == 1

      # Not so simple arithmetic: applying starting on the left:
      # ϱ₁₂*ϱ₂₁⁻¹*λ₁₂*ε₂ == σ₂₁₃₄

      g = A(Groups.transvection_R(1,2))
      x1, x2, x3, x4 = Groups.domain(A)
      @test g(Groups.domain(A)) == (x1*x2, x2, x3, x4)
      g = g*inv(A(Groups.transvection_R(2,1)))
      @test g(Groups.domain(A)) == (x1*x2, x1^-1, x3, x4)
      g = g*A(Groups.transvection_L(1,2))
      @test g(Groups.domain(A)) == (x2, x1^-1, x3, x4)
      g = g*A(Groups.flip(2))
      @test g(Groups.domain(A)) == (x2, x1, x3, x4)

      @test g(Groups.domain(A)) == A(Groups.AutSymbol(perm"(1,2)(4)"))(Groups.domain(A))

      @test g == A(Groups.AutSymbol(perm"(1,2)(4)"))

      g_im = g(Groups.domain(A))
      @test length.(g_im) == (1,1,1,1)
   end

   @testset "specific Aut(F4) tests" begin
      N = 4
      G = AutGroup(FreeGroup(N))
      S = G.gens
      @test isa(S, Vector{Groups.AutSymbol})
      S = [G(s) for s in unique(S)]
      @test isa(S, Vector{Automorphism{N}})
      @test S == gens(G)
      @test length(S) == 51
      S_inv = [S..., [inv(s) for s in S]...]
      @test length(unique(S_inv)) == 75

      G = AutGroup(FreeGroup(N), special=true)
      S = gens(G)
      S_inv = [one(G), S..., [inv(s) for s in S]...]
      S_inv = unique(S_inv)
      B_2 = [i*j for (i,j) in Base.product(S_inv, S_inv)]
      @test length(B_2) == 2401
      @test length(unique(B_2)) == 1777
   end

   @testset "linear_repr tests" begin
      N = 3
      G = AutGroup(FreeGroup(N))
      S = unique([gens(G); inv.(gens(G))])
      R = 3

      @test Groups.linear_repr(one(G)) isa Matrix{Int}
      @test Groups.linear_repr(one(G)) == Matrix{Int}(I, N, N)

      M = Matrix{Int}(I, N, N)
      M[1,2] = 1
      ϱ₁₂ = G(Groups.transvection_R(1,2))
      λ₁₂ = G(Groups.transvection_R(1,2))

      @test Groups.linear_repr(ϱ₁₂) == M
      @test Groups.linear_repr(λ₁₂) == M

      M[1,2] = -1

      @test Groups.linear_repr(ϱ₁₂^-1) == M
      @test Groups.linear_repr(λ₁₂^-1) == M

      @test Groups.linear_repr(ϱ₁₂*λ₁₂^-1) == Matrix{Int}(I, N, N)
      @test Groups.linear_repr(λ₁₂^-1*ϱ₁₂) == Matrix{Int}(I, N, N)

      M = Matrix{Int}(I, N, N)
      M[2,2] = -1
      ε₂ = G(Groups.flip(2))

      @test Groups.linear_repr(ε₂) == M
      @test Groups.linear_repr(ε₂^2) == Matrix{Int}(I, N, N)

      M = [0 1 0; 0 0 1; 1 0 0]

      σ = G(Groups.AutSymbol(perm"(1,2,3)"))
      @test Groups.linear_repr(σ) == M
      @test Groups.linear_repr(σ^3) == Matrix{Int}(I, 3, 3)
      @test Groups.linear_repr(σ)^3 == Matrix{Int}(I, 3, 3)

      function test_homomorphism(S, r)
         for elts in Iterators.product([[g for g in S] for _ in 1:r]...)
            prod(Groups.linear_repr.(elts)) == Groups.linear_repr(prod(elts)) || error("linear representaton test failed at $elts")
         end
         return 0
      end

      @test test_homomorphism(S, R) == 0
   end
end
