@testset "Automorphisms" begin
   using Nemo
   G = PermutationGroup(4)

   @testset "AutSymbol" begin
      @test_throws MethodError Groups.AutSymbol("a")
      @test_throws MethodError Groups.AutSymbol("a", 1)
      f = AutSymbol("a", 1, :(a()), v -> v)
      @test isa(f, Groups.GSymbol)
      @test isa(f, Groups.AutSymbol)
      @test isa(Groups.perm_autsymbol(G([1,2,3,4])), Groups.AutSymbol)
      @test isa(Groups.rmul_autsymbol(1,2), Groups.AutSymbol)
      @test isa(Groups.lmul_autsymbol(3,4), Groups.AutSymbol)
      @test isa(Groups.flip_autsymbol(3), Groups.AutSymbol)
   end

   a,b,c,d = generators(FreeGroup(4))
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
      f = AutSymbol("a", 1, :(a()), v -> v)
      @test isa(GWord(f), GWord)
      @test isa(GWord(f), AutGroupElem)
      @test isa(AutGroupElem(f), AutGroupElem)
      @test isa(AutGroup(FreeGroup(3)), Group)
      @test isa(AutGroup(FreeGroup(1)), FPGroup)
      A = AutGroup(FreeGroup(1))
      @test isa(f*f, AutWord)
      @test isa(f^2, AutWord)
      @test isa(f^-1, AutWord)
   end
#
#    @testset "eltary functions" begin
#       f = perm_autsymbol([2,1,4,3])
#       @test isa(inv(f), AutSymbol)
#       @test isa(f^-1, AutWord)
#       @test f^-1 == GWord(inv(f))
#       @test inv(f) == f
#    end
#
#    @testset "reductions/arithmetic" begin
#       f = perm_autsymbol([2,1,4,3])
#       f² = Groups.r_multiply(AutWord(f), [f], reduced=false)
#       @test Groups.simplify_perms!(f²) == false
#       @test f² == one(typeof(f*f))
#
#       a = rmul_autsymbol(1,2)*flip_autsymbol(2)
#       b = flip_autsymbol(2)*inv(rmul_autsymbol(1,2))
#       @test a*b == b*a
#       @test a^3 * b^3 == one(a)
#    end
#
#    @testset "specific Aut(𝔽₄) tests" begin
#       N = 4
#       import Combinatorics.nthperm
#       SymmetricGroup(n) = [nthperm(collect(1:n), k) for k in 1:factorial(n)]
#       indexing = [[i,j] for i in 1:N for j in 1:N if i≠j]
#
#       σs = [perm_autsymbol(perm) for perm in SymmetricGroup(N)[2:end]];
#       ϱs = [rmul_autsymbol(i,j) for (i,j) in indexing]
#       λs = [lmul_autsymbol(i,j) for (i,j) in indexing]
#       ɛs = [flip_autsymbol(i) for i in 1:N];
#
#       S = vcat(ϱs, λs, σs, ɛs)
#       S = vcat(S, [inv(s) for s in S])
#       @test isa(S, Vector{AutSymbol})
#       @test length(S) == 102
#       @test length(unique(S)) == 75
#       S₁ = [GWord(s) for s in unique(S)]
#       @test isa(S₁, Vector{AutWord})
#       p = prod(S₁)
#       @test length(p) == 53
#    end

end
