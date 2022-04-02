using BenchmarkTools
using Test

using Groups

function wl_ball(F; radius::Integer)
    g, state = iterate(F)
    while length(word(g)) <= radius
        res = iterate(F, state)
        isnothing(res) && break
        g, state = res
    end
    elts = collect(state.seen)
    elts = resize!(elts, length(elts)-1)
    return elts
end

@testset "Benchmarks" begin
    N = 4

    @testset "iteration: FreeGroup" begin
        FN = FreeGroup(N)
        R = 8

        let G = FN
            S = unique([gens(G); inv.(gens(G))])

            sizes1 = last(Groups.wlmetric_ball(S, radius=R, threading=false))
            sizes2 = last(Groups.wlmetric_ball(S, radius=R, threading=true))

            l = length(wl_ball(G, radius=R))

            @test sizes1 == sizes2
            @test last(sizes1) == l

            @info "Ball of radius $R in $(parent(first(S)))" sizes=sizes1
            @info "serial"
            @time Groups.wlmetric_ball(S, radius=R, threading=false)
            @info "threaded"
            @time Groups.wlmetric_ball(S, radius=R, threading=true)
            @info "iteration"
            @time wl_ball(G, radius=R)
        end
    end

    @testset "iteration: SAut(F_n)" begin
        R = 4
        FN = FreeGroup(N)
        SAutFN = SpecialAutomorphismGroup(FN)

        let G = SAutFN
            S = unique([gens(G); inv.(gens(G))])

            sizes1 = last(Groups.wlmetric_ball(S, radius=R, threading=false))
            sizes2 = last(Groups.wlmetric_ball(S, radius=R, threading=true))

            l = length(wl_ball(G, radius=R))

            @test sizes1 == sizes2
            @test last(sizes1) == l

            @info "Ball of radius $R in $(parent(first(S)))" sizes=sizes1
            @info "serial"
            @time Groups.wlmetric_ball(S, radius=R, threading=false)
            @info "threaded"
            @time Groups.wlmetric_ball(S, radius=R, threading=true)
            @info "iteration"
            @time wl_ball(G, radius=R)
        end
    end
end
