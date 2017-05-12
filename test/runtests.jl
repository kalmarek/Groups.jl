using Groups
using Base.Test

# write your own tests here

@testset "Groups" begin

   @testset "Groups.FreeSymbols" begin
      s = Groups.FreeSymbol("s")
      t = Groups.FreeSymbol("t")
      @testset "defines" begin
         @test isa(Groups.FreeSymbol("aaaaaaaaaaaaaaaa"), Groups.GSymbol)
         @test Groups.FreeSymbol("abc").pow == 1
         @test isa(s, Groups.FreeSymbol)
         @test isa(t, Groups.FreeSymbol)
      end
      @testset "eltary functions" begin
         @test length(s) == 1
         @test Groups.change_pow(s, 0) == Groups.change_pow(t, 0)
         @test length(Groups.change_pow(s, 0)) == 0
         @test inv(s).pow == -1
         @test Groups.FreeSymbol("s", 3) == Groups.change_pow(s, 3)
         @test Groups.FreeSymbol("s", 3) != Groups.FreeSymbol("t", 3)
         @test Groups.change_pow(inv(s), -3) == inv(Groups.change_pow(s, 3))
      end
      @testset "powers" begin
         s⁴ = Groups.change_pow(s,4)
         @test s⁴.pow == 4
         @test Groups.change_pow(s, 4) == Groups.FreeSymbol("s", 4)
      end
   end

   @testset "FreeGroupElems" begin
      s = Groups.FreeSymbol("s")
      t = Groups.FreeSymbol("t", -2)
      @testset "defines" begin
         @test isa(Groups.GWord(s), Groups.GWord)
         @test isa(Groups.GWord(s), FreeGroupElem)
         @test isa(FreeGroupElem(s), Groups.GWord)
         @test isa(convert(FreeGroupElem, s), Groups.GWord)
         @test isa(convert(FreeGroupElem, s), FreeGroupElem)
         @test isa(Vector{FreeGroupElem}([s,t]), Vector{FreeGroupElem})
         @test length(FreeGroupElem(s)) == 1
         @test length(FreeGroupElem(t)) == 2
      end

      @testset "eltary functions" begin
         @test_skip (s*s).symbols == (s^2).symbols
         @test_skip Vector{Groups.GWord{Groups.FreeSymbol}}([s,t]) ==
            Vector{FreeGroupElem}([s,t])
         @test_skip Vector{Groups.GWord}([s,t]) ==
            [Groups.GWord(s), Groups.GWord(t)]
         @test_skip hash([t^1,s^1]) == hash([t^2*inv(t),s*inv(s)*s])
      end
   end

   @testset "FreeGroup" begin
      @test isa(FreeGroup(["s", "t"]), Nemo.Group)
      G = FreeGroup(["s", "t"])

      @testset "elements constructors" begin
         @test isa(G(), FreeGroupElem)
         @test eltype(G.gens) == Groups.FreeSymbol
         @test length(G.gens) == 2
         @test_skip eltype(G.rels) == FreeGroupElem
         @test_skip length(G.rels) == 0
         @test eltype(generators(G)) == FreeGroupElem
         @test length(generators(G)) == 2
      end

      s, t = generators(G)

      @testset "internal arithmetic" begin
         t_symb = Groups.FreeSymbol("t")
         tt = deepcopy(t)
         @test string(Groups.r_multiply!(tt,[inv(t_symb)]; reduced=true)) ==
            "(id)"
         tt = deepcopy(t)
         @test string(Groups.r_multiply!(tt,[inv(t_symb)]; reduced=false)) ==
            "t*t^-1"
         tt = deepcopy(t)
         @test string(Groups.l_multiply!(tt,[inv(t_symb)]; reduced=true)) ==
            "(id)"
         tt = deepcopy(t)
         @test string(Groups.l_multiply!(tt,[inv(t_symb)]; reduced=false)) ==
            "t^-1*t"
      end

      @testset "reductions" begin
         @test length(G().symbols) == 1
         @test length((G()*G()).symbols) == 0
         @test G() == G()*G()
         w = deepcopy(s)
         push!(w.symbols, (s^-1).symbols[1])
         @test Groups.reduce!(w) == parent(w)()
         o = (t*s)^3
         @test o == t*s*t*s*t*s
         p = (t*s)^-3
         @test p == s^-1*t^-1*s^-1*t^-1*s^-1*t^-1
         @test o*p == parent(o*p)()
         w = FreeGroupElem([o.symbols..., p.symbols...])
         w.parent = G
         @test Groups.reduce!(w).symbols ==Vector{Groups.FreeSymbol}([])
      end

      @testset "binary/inv operations" begin
         @test parent(s) == G
         @test parent(s) === parent(deepcopy(s))
         @test isa(s*t, FreeGroupElem)
         @test parent(s*t) == parent(s^2)
         @test s*s == s^2
         @test inv(s*s) == inv(s^2)
         @test inv(s)^2 == inv(s^2)
         @test inv(s)*inv(s) == inv(s^2)
         @test inv(s*t) == inv(t)*inv(s)
         w = s*t*s^-1
         @test inv(w) == s*t^-1*s^-1
         @test (t*s*t^-1)^10 == t*s^10*t^-1
         @test (t*s*t^-1)^-10 == t*s^-10*t^-1
      end

      @testset "replacements" begin
         a = Groups.FreeSymbol("a")
         b = Groups.FreeSymbol("b")
         @test Groups.is_subsymbol(a, Groups.change_pow(a,2)) == true
         @test Groups.is_subsymbol(a, Groups.change_pow(a,-2)) == false
         @test Groups.is_subsymbol(b, Groups.change_pow(a,-2)) == false
         @test Groups.is_subsymbol(inv(b), Groups.change_pow(b,-2)) == true
         c = s*t*s^-1*t^-1
         @test findfirst(c, s^-1*t^-1) == 3
         @test findnext(c*s^-1, s^-1*t^-1,3) == 3
         @test findnext(c*s^-1*t^-1, s^-1*t^-1,4) == 5
         @test findfirst(c*t, c) == 0
         w = s*t*s^-1
         subst = Dict{FreeGroupElem, FreeGroupElem}(w => s^1, s*t^-1 => t^4)
         @test Groups.replace(c, 1, s*t, G()) == s^-1*t^-1
         @test Groups.replace(c, 1, w, subst[w]) == s*t^-1
         @test Groups.replace(s*c*t^-1, 1, w, subst[w]) == s^2*t^-2
         @test Groups.replace(t*c*t, 2, w, subst[w]) == t*s
         @test Groups.replace_all!(s*c*s*c*s, subst) == s*t^4*s*t^4*s
      end
   end


   @testset "Automorphisms" begin
      @testset "AutSymbol" begin
         @test_throws MethodError AutSymbol("a")
         @test_throws MethodError AutSymbol("a", 1)
         f = AutSymbol("a", 1, :(a()), v -> v)
         @test isa(f, Groups.GSymbol)
         @test isa(f, Groups.AutSymbol)
         @test isa(symmetric_AutSymbol([1,2,3,4]), AutSymbol)
         @test isa(rmul_AutSymbol(1,2), AutSymbol)
         @test isa(lmul_AutSymbol(3,4), AutSymbol)
         @test isa(flip_AutSymbol(3), AutSymbol)
      end

   #    @testset "flip_AutSymbol correctness" begin
   #       a,b,c,d = [FreeGroupElem(Groups.FreeSymbol(i)) for i in ["a", "b", "c", "d"]]
   #       domain = [a,b,c,d]
   #       @test flip_AutSymbol(1)(domain) == [a^-1, b,c,d]
   #       @test flip_AutSymbol(2)(domain) == [a, b^-1,c,d]
   #       @test flip_AutSymbol(3)(domain) == [a, b,c^-1,d]
   #       @test flip_AutSymbol(4)(domain) == [a, b,c,d^-1]
   #       @test inv(flip_AutSymbol(1))(domain) == [a^-1, b,c,d]
   #       @test inv(flip_AutSymbol(2))(domain) == [a, b^-1,c,d]
   #       @test inv(flip_AutSymbol(3))(domain) == [a, b,c^-1,d]
   #       @test inv(flip_AutSymbol(4))(domain) == [a, b,c,d^-1]
   #    end
   #
   #    @testset "symmetric_AutSymbol correctness" begin
   #       a,b,c,d = [FreeGroupElem(Groups.FreeSymbol(i)) for i in ["a", "b", "c", "d"]]
   #       domain = [a,b,c,d]
   #       σ = symmetric_AutSymbol([1,2,3,4])
   #       @test σ(domain) == domain
   #       @test inv(σ)(domain) == domain
   #
   #       σ = symmetric_AutSymbol([2,3,4,1])
   #       @test σ(domain) == [b, c, d, a]
   #       @test inv(σ)(domain) == [d, a, b, c]
   #
   #       σ = symmetric_AutSymbol([2,1,4,3])
   #       @test σ(domain) == [b, a, d, c]
   #       @test inv(σ)(domain) == [b, a, d, c]
   #
   #       σ = symmetric_AutSymbol([2,3,1,4])
   #       @test σ(domain) == [b,c,a,d]
   #       @test inv(σ)(domain) == [c,a,b,d]
   #    end
   #
   #    @testset "mul_AutSymbol correctness" begin
   #       a,b,c,d = [FreeGroupElem(Groups.FreeSymbol(i)) for i in ["a", "b", "c", "d"]]
   #       domain = [a,b,c,d]
   #       i,j = 1,2
   #       r = rmul_AutSymbol(i,j)
   #       l = lmul_AutSymbol(i,j)
   #       @test r(domain) == [a*b,b,c,d]
   #       @test inv(r)(domain) == [a*b^-1,b,c,d]
   #       @test l(domain) == [b*a,b,c,d]
   #       @test inv(l)(domain) == [b^-1*a,b,c,d]
   #
   #       i,j = 3,1
   #       r = rmul_AutSymbol(i,j)
   #       l = lmul_AutSymbol(i,j)
   #       @test r(domain) == [a,b,c*a,d]
   #       @test inv(r)(domain) == [a,b,c*a^-1,d]
   #       @test l(domain) == [a,b,a*c,d]
   #       @test inv(l)(domain) == [a,b,a^-1*c,d]
   #
   #
   #       i,j = 4,3
   #       r = rmul_AutSymbol(i,j)
   #       l = lmul_AutSymbol(i,j)
   #       @test r(domain) == [a,b,c,d*c]
   #       @test inv(r)(domain) == [a,b,c,d*c^-1]
   #       @test l(domain) == [a,b,c,c*d]
   #       @test inv(l)(domain) == [a,b,c,c^-1*d]
   #
   #
   #       i,j = 2,4
   #       r = rmul_AutSymbol(i,j)
   #       l = lmul_AutSymbol(i,j)
   #       @test r(domain) == [a,b*d,c,d]
   #       @test inv(r)(domain) == [a,b*d^-1,c,d]
   #       @test l(domain) == [a,d*b,c,d]
   #       @test inv(l)(domain) == [a,d^-1*b,c,d]
   #    end
   #
   #    @testset "AutWords" begin
   #       f = AutSymbol("a", 1, :(a()), v -> v)
   #       @test isa(GWord(f), GWord)
   #       @test isa(GWord(f), AutWord)
   #       @test isa(AutWord(f), AutWord)
   #       @test isa(f*f, AutWord)
   #       @test isa(f^2, AutWord)
   #       @test isa(f^-1, AutWord)
   #    end
   #
   #    @testset "eltary functions" begin
   #       f = symmetric_AutSymbol([2,1,4,3])
   #       @test isa(inv(f), AutSymbol)
   #       @test isa(f^-1, AutWord)
   #       @test f^-1 == GWord(inv(f))
   #       @test inv(f) == f
   #    end
   #
   #    @testset "reductions/arithmetic" begin
   #       f = symmetric_AutSymbol([2,1,4,3])
   #       f² = Groups.r_multiply(AutWord(f), [f], reduced=false)
   #       @test Groups.simplify_perms!(f²) == false
   #       @test f² == one(typeof(f*f))
   #
   #       a = rmul_AutSymbol(1,2)*flip_AutSymbol(2)
   #       b = flip_AutSymbol(2)*inv(rmul_AutSymbol(1,2))
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
   #       σs = [symmetric_AutSymbol(perm) for perm in SymmetricGroup(N)[2:end]];
   #       ϱs = [rmul_AutSymbol(i,j) for (i,j) in indexing]
   #       λs = [lmul_AutSymbol(i,j) for (i,j) in indexing]
   #       ɛs = [flip_AutSymbol(i) for i in 1:N];
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
   # end

end
