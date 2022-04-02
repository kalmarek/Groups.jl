using Test
using Groups
using PermutationGroups

import Logging

import KnuthBendix: Word

using GroupsCore
include(joinpath(pathof(GroupsCore), "..", "..", "test", "conformance_test.jl"))

@testset "Groups" begin

    include("free_groups.jl")
    include("fp_groups.jl")

    include("matrix_groups.jl")
    include("AutFn.jl")

    include("homomorphisms.jl")

    include("AutSigma_41.jl")
    include("AutSigma3.jl")

    include("group_constructions.jl")

    # if !haskey(ENV, "CI")
    #    include("benchmarks.jl")
    # end
end
