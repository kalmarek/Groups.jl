using BenchmarkTools
using Test

using Groups

function wl_ball(F; radius::Integer)
    g, state = iterate(F)
    sizes = Int[]
    while length(sizes) ≤ radius
        res = iterate(F, state)
        isnothing(res) && break
        g, state = res
        if length(word(g)) > length(sizes)
            push!(sizes, length(state.seen) - 1)
        end
    end
    elts = collect(state.seen)
    resize!(elts, sizes[end] - 1)
    return elts, sizes[2:end]
end

@testset "Benchmarks" begin
    N = 4

    @testset "iteration: FreeGroup" begin
        FN = FreeGroup(N)
        R = 8

        let G = FN
            S = unique([gens(G); inv.(gens(G))])

            sizes1 = last(Groups.wlmetric_ball(S; radius = R))
            sizes2 = last(wl_ball(G; radius = R))

            @test sizes1 == sizes2

            @info "Ball of radius $R in $(parent(first(S)))" sizes = sizes1
            @info "serial"
            @time Groups.wlmetric_ball(S, radius = R)
            @info "iteration"
            @time wl_ball(G, radius = R)
        end
    end

    @testset "iteration: SAut(F_n)" begin
        R = 4
        FN = FreeGroup(N)
        SAutFN = SpecialAutomorphismGroup(FN)

        let G = SAutFN
            S = unique([gens(G); inv.(gens(G))])

            sizes1 = last(Groups.wlmetric_ball(S; radius = R))
            sizes2 = last(wl_ball(G; radius = R))

            @test sizes1 == sizes2

            @info "Ball of radius $R in $(parent(first(S)))" sizes = sizes1
            @info "serial"
            @time Groups.wlmetric_ball(S, radius = R)
            @info "iteration"
            @time wl_ball(G, radius = R)
        end
    end
end
