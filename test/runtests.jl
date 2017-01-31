using Groups
using Base.Test

# write your own tests here
s = FGSymbol("s")
t = FGSymbol("t")

@testset "FGSymbols" begin
    @testset "defines" begin
        @test isa(FGSymbol(string(Char(rand(50:2000)))), Groups.GSymbol)
        @test FGSymbol("abc").pow == 1
        @test isa(s, FGSymbol)
        @test isa(t, FGSymbol)
    end
    @testset "eltary functions" begin
        @test length(s) == 1
        @test one(s) == s^0
        @test one(s) == one(FGSymbol)
        @test Groups.change_pow(s,0) == one(s)
        @test length(one(s)) == 0
        @test inv(s).pow == -1
        @test FGSymbol("s", 3) == Groups.change_pow(s,3)
        @test s^2 ≠ t^2

    end
    @testset "powers" begin
        s⁴ = Groups.change_pow(s,4)
        @test s⁴.pow == 4
        @test (s^4).symbols[1] == Groups.change_pow(s,4)
        @test s*s == s^2
        @test inv(s*s) == inv(s^2)
        @test inv(s)^2 == inv(s^2)
        @test inv(s)*inv(s) == inv(s^2)
        @test inv(s*s) == inv(s)*inv(s)
    end
end

@testset "FGWords" begin
    @testset "defines" begin
        @test isa(Groups.GWord(s), Groups.GWord)
        @test isa(Groups.GWord(s), FGWord)
        @test isa(FGWord(s), Groups.GWord)
        @test isa(convert(FGWord, s), GWord)
        @test isa(convert(FGWord, s), FGWord)
        @test isa(Vector{FGWord}([s,t]), Vector{FGWord})
        @test Vector{GWord{FGSymbol}}([s,t]) == Vector{FGWord}([s,t])
        @test isa(s*s, FGWord)
        @test s*s == s^2
        @test t*s ≠ s*t
        @test Vector{GWord}([s,t]) == [s^2*s^-1, t]
        @test hash([t^1,s^1]) == hash([t^2*inv(t),s*inv(s)*s])
    end
    @testset "eltary functions" begin
        @test length(FGWord(s)) == 1
        @test length(s*s) == 2
        @test length(s*s^-1) == 0
        @test length(s*t^-1) == 2
        @test isa(one(FGWord), FGWord)
        @test one(FGWord).symbols == Vector{FGSymbol}([one(FGSymbol)])
        @test isa(one(Groups.GWord{FGSymbol}), Groups.GWord{FGSymbol})
        w = s*t*s^-1
        @test isa(one(w), FGWord)
        @test inv(s*t) == t^-1*s^-1
        @test inv(w) == s*t^-1*s^-1

    end

    @testset "reductions" begin
        @test one(FGWord) == one(s)*one(s)
        w = GWord{FGSymbol}([s])
        push!(w.symbols, (s^-1).symbols[1])
        @test Groups.reduce!(w) == one(FGWord)
        o = (t*s)^3
        @test o == t*s*t*s*t*s
        p = (t*s)^-3
        @test p == s^-1*t^-1*s^-1*t^-1*s^-1*t^-1
        @test o*p == one(FGWord)
        w = FGWord([o.symbols..., p.symbols...])
        @test Groups.reduce!(w).symbols ==Vector{FGSymbol}([])
    end
    @testset "arithmetic" begin
        @test Groups.r_multiply!(FGWord(t),[s,t]; reduced=true) == t*s*t
        @test Groups.r_multiply!(FGWord(t),[s,t]; reduced=false) == t*s*t

        @test Groups.l_multiply!(FGWord(t),[s,t]; reduced=true) == t*s*t
        @test Groups.l_multiply!(FGWord(t),[s,t]; reduced=false) == t*s*t
        @test (t*s*t^-1)^10 == t*s^10*t^-1
        @test (t*s*t^-1)^-10 == t*s^-10*t^-1
    end

    @testset "replacements" begin
        @test Groups.is_subsymbol(s, Groups.change_pow(s,2)) == true
        @test Groups.is_subsymbol(s, Groups.change_pow(s,-2)) == false
        @test Groups.is_subsymbol(t, Groups.change_pow(s,-2)) == false
        @test Groups.is_subsymbol(inv(t), Groups.change_pow(t,-2)) == true
        c = s*t*s^-1*t^-1
        @test findfirst(c, s^-1*t^-1) == 3
        @test findnext(c*s^-1, s^-1*t^-1,3) == 3
        @test findnext(c*s^-1*t^-1, s^-1*t^-1,4) == 5
        @test findfirst(c*t, c) == 0
        w = s*t*s^-1
        subst = Dict{FGWord, FGWord}(w => s^1, s*t^-1 => t^4)
        @test Groups.replace(c, 1, s*t, one(FGWord)) == s^-1*t^-1

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
        f = AutSymbol("a", 1, :(a(0)), v -> v, v -> v)
        @test isa(f, GSymbol)
        @test isa(f, AutSymbol)
        @test isa(symmetric_AutSymbol([1,2,3,4]), AutSymbol)
        @test isa(rmul_AutSymbol(1,2), AutSymbol)
        @test isa(lmul_AutSymbol(3,4), AutSymbol)
        @test isa(flip_AutSymbol(3), AutSymbol)
    end

    @testset "AutWords" begin
        f = AutSymbol("a", 1, :(a(0)), v -> v, v -> v)
        @test isa(GWord(f), GWord)
        @test isa(GWord(f), AutWord)
        @test isa(AutWord(f), AutWord)
        @test isa(f*f, AutWord)
        @test isa(f^2, AutWord)
        @test isa(f^-1, AutWord)
    end
    @testset "eltary functions" begin
        f = symmetric_AutSymbol([2,1,4,3])
        @test isa(inv(f), AutSymbol)
        @test isa(f^-1, AutWord)
        @test f^-1 == GWord(inv(f))
        @test inv(f) == f
    end
    @testset "reductions/arithmetic" begin
        f = symmetric_AutSymbol([2,1,4,3])
        f² = Groups.r_multiply(AutWord(f), [f], reduced=false)
        @test Groups.simplify_perms!(f²) == false
        @test f² == one(typeof(f*f))

        a = rmul_AutSymbol(1,2)*flip_AutSymbol(2)
        b = flip_AutSymbol(2)*inv(rmul_AutSymbol(1,2))
        @test a*b == b*a
        @test a^3 * b^3 == one(a)
    end
    @testset "specific Aut(𝔽₄) tests" begin
        N = 4
        import Combinatorics.nthperm
        SymmetricGroup(n) = [nthperm(collect(1:n), k) for k in 1:factorial(n)]
        indexing = [[i,j] for i in 1:N for j in 1:N if i≠j]

        σs = [symmetric_AutSymbol(perm) for perm in SymmetricGroup(N)[2:end]];
        ϱs = [rmul_AutSymbol(i,j) for (i,j) in indexing]
        λs = [lmul_AutSymbol(i,j) for (i,j) in indexing]
        ɛs = [flip_AutSymbol(i) for i in 1:N];

        S = vcat(ϱs, λs, σs, ɛs)
        S = vcat(S, [inv(s) for s in S])
        @test isa(S, Vector{AutSymbol})
        @test length(S) == 102
        @test length(unique(S)) == 75
        S₁ = [GWord(s) for s in unique(S)]
        @test isa(S₁, Vector{AutWord})
        p = prod(S₁)
        @test length(p) == 75
        @test Groups.simplify_perms!(p) == false
        @test length(p) == 53
        @test Groups.join_free_symbols!(p) == true


    end
end
