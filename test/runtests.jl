using Test
using Groups
using PermutationGroups

import Logging

import KnuthBendix: Word

using GroupsCore
include(joinpath(pathof(GroupsCore), "..", "..", "test", "conformance_test.jl"))

@testset "Groups" begin

    _, t = @timed include("free_groups.jl")
    @info "free_groups.jl took $(round(t, digits=2))s"
    _, t = @timed include("fp_groups.jl")
    @info "fp_groups.jl took $(round(t, digits=2))s"

    _, t = @timed include("matrix_groups.jl")
    @info "matrix_groups.jl took $(round(t, digits=2))s"
    _, t = @timed include("AutFn.jl")
    @info "AutFn.jl took $(round(t, digits=2))s"

    _, t = @timed include("homomorphisms.jl")
    @info "homomorphisms.jl took $(round(t, digits=2))s"

    if !haskey(ENV, "CI")
        _, t = @timed include("AutSigma_41.jl")
        @info "AutSigma_41 took $(round(t, digits=2))s"
        _, t = @timed include("AutSigma3.jl")
        @info "AutSigma3 took $(round(t, digits=2))s"
    end

    _, t = @timed include("group_constructions.jl")
    @info "Constructions took $(round(t, digits=2))s"
end

if !haskey(ENV, "CI")
    include("benchmarks.jl")
end
