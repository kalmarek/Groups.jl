using Test
using AbstractAlgebra
using Groups

using LinearAlgebra

@testset "Groups" begin

   @testset "generate balls" begin
      M = MatrixAlgebra(zz, 3)
      w = one(M); w[1,2] = 1;
      r = one(M); r[2,3] = -3;
      s = one(M); s[1,3] = 2; s[3,2] = -1;

      S = [w,r,s]; S = unique([S; inv.(S)]);
      _, sizes = Groups.generate_balls(S, radius=4);
      @test sizes == [7, 33, 141, 561]
   end

   include("FreeGroup-tests.jl")
   include("AutGroup-tests.jl")
   include("DirectPower-tests.jl")
   include("WreathProd-tests.jl")
   include("FPGroup-tests.jl")
end
