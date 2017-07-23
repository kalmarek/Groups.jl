using Groups
using Base.Test

# write your own tests here

@testset "Groups" begin
   include("FreeGroup-tests.jl")
   include("AutGroup-tests.jl")
   include("DirectProd-tests.jl")
end
