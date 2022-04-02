using Groups.MatrixGroups

@testset "Matrix Groups" begin

    @testset "SL(n, ℤ)" begin
        SL3Z = SpecialLinearGroup{3}(Int8)

        S = gens(SL3Z); union!(S, inv.(S))

        E, sizes = Groups.wlmetric_ball(S, radius=4)

        @test sizes = [13, 121, 883, 5455]

        @testset "GroupsCore conformance" begin
            test_Group_interface(SL3Z)
            g = A(rand(1:length(alphabet(SL3Z)), 10))
            h = A(rand(1:length(alphabet(SL3Z)), 10))

            test_GroupElement_interface(g, h)
        end

    end

end
