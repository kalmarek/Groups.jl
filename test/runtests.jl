using Test
using AbstractAlgebra
using Groups

using LinearAlgebra

@testset "Groups" begin
   include("FreeGroup-tests.jl")
   include("AutGroup-tests.jl")
   include("DirectProd-tests.jl")
   include("WreathProd-tests.jl")
end
