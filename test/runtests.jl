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

@testset "GWords" begin
    @testset "defines" begin
        @test isa(Groups.GWord(s), Groups.GWord)
        @test isa(Groups.GWord(s), FGWord)
        @test isa(FGWord(s), Groups.GWord)
        @test isa(s*s, FGWord)
        @test s*s == s^2
        @test t*s ≠ s*t
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
        @test Groups.freegroup_reduce!(w) == one(FGWord)
        o = (t*s)^3
        @test o == t*s*t*s*t*s
        p = (t*s)^-3
        @test p == s^-1*t^-1*s^-1*t^-1*s^-1*t^-1
        @test o*p == one(FGWord)
        w = FGWord([o.symbols..., p.symbols...])
        @test Groups.freegroup_reduce!(w).symbols ==Vector{FGSymbol}([])
    end
    @testset "arithmetic" begin
        @test Groups.r_multiply!(FGWord(t),[s,t]; reduced=true) == t*s*t
        @test Groups.r_multiply!(FGWord(t),[s,t]; reduced=false) == t*s*t

        @test Groups.l_multiply!(FGWord(t),[s,t]; reduced=true) == t*s*t
        @test Groups.l_multiply!(FGWord(t),[s,t]; reduced=false) == t*s*t
        @test (t*s*t^-1)^10 == t*s^10*t^-1
        @test (t*s*t^-1)^-10 == t*s^-10*t^-1
    end
@testset "Automorphisms" begin
    @testset "AutSymbol" begin
        @test_throws MethodError AutSymbol("a")
        @test_throws MethodError AutSymbol("a", 1)
        f = AutSymbol("a", 1, :(a(0)))
        @test isa(f, GSymbol)
        @test isa(f, AutSymbol)
        @test isa(symmetric_AutSymbol([1,2,3,4]), AutSymbol)
        @test isa(rmul_AutSymbol(1,2), AutSymbol)
        @test isa(lmul_AutSymbol(3,4), AutSymbol)
        @test isa(flip_AutSymbol(3), AutSymbol)
    end

    @testset "AutWords" begin
        f = AutSymbol("a", 1, :(a(0)))
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
end
