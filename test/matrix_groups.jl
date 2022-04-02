using Groups.MatrixGroups

@testset "Matrix Groups" begin
    @testset "SL(n, ℤ)" begin
        SL3Z = SpecialLinearGroup{3}(Int8)

        S = gens(SL3Z); union!(S, inv.(S))

        E, sizes = Groups.wlmetric_ball(S, radius=4)

        @test sizes == [13, 121, 883, 5455]

        E(i,j) = SL3Z([A[MatrixGroups.ElementaryMatrix{3}(i,j, Int8(1))]])

        A = alphabet(SL3Z)
        w = E(1,2)
        r = E(2,3)^-3
        s = E(1,3)^2*E(3,2)^-1

        S = [w,r,s]; S = unique([S; inv.(S)]);
        _, sizes = Groups.wlmetric_ball(S, radius=4);
        @test sizes == [7, 33, 141, 561]
        _, sizes = Groups.wlmetric_ball_serial(S, radius=4);
        @test sizes == [7, 33, 141, 561]

        @testset "GroupsCore conformance" begin
            test_Group_interface(SL3Z)
            g = A(rand(1:length(alphabet(SL3Z)), 10))
            h = A(rand(1:length(alphabet(SL3Z)), 10))

            test_GroupElement_interface(g, h)
        end
    end
end
