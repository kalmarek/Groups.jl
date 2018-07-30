using Base.Test
using AbstractAlgebra
using Groups

@testset "Groups" begin
   include("FreeGroup-tests.jl")
   include("AutGroup-tests.jl")
   include("DirectProd-tests.jl")
   include("WreathProd-tests.jl")
end
