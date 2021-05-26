@testset "New.FreeGroup" begin

    A3 = Alphabet([:a, :b, :c, :A, :B, :C], [4,5,6,1,2,3])
    F3 = New.FreeGroup([:a, :b, :c], A3)
    @test F3 isa New.FreeGroup

    @test gens(F3) isa Vector

    @test eltype(F3) <: New.FPGroupElement{<:New.FreeGroup}

    w = F3([1,2,3,4])
    W = inv(w)
    @test deepcopy(w) !== w
    @test deepcopy(w).word !== w.word

    @test isone(w*W)

    @test New.alphabet(w) == A3

    @testset "generic iteration" begin
        w, s = iterate(F3)
        @test isone(w)
        w, s = iterate(F3, s)
        @test w == gens(F3, 1)

        a,b,c = gens(F3)

        function test_iteration(F, n=1000)
            w, s = iterate(F)
            for i in 1:n
                w, s = iterate(F, s)
            end
            return w
        end

        k = test_iteration(F3, 10)
        @test k == a*b^-1

        @time k = test_iteration(F3, 1000)
        @test k == (a^2)*c^2*a^-1
    end

    @testset "wl_ball" begin
        function wl_ball(F; radius::Integer)
            g, state = iterate(F)
            while length(New.word(g)) <= radius
                res = iterate(F, state)
                isnothing(res) && break
                g, state = res
            end
            elts = collect(state.seen)
            elts = resize!(elts, length(elts)-1)
            return elts
        end

        E4 = wl_ball(F3, radius=4)
        @test length(E4) == 937
        @test New.word(last(E4)) == Word([6])^4

        E8, t, _ = @timed wl_ball(F3, radius=8)
        @test length(E8) == 585937
        @test New.word(last(E8)) == Word([6])^8
        @test t/10^9 < 1
    end

    @testset "GroupsCore conformance" begin
        test_Group_interface(F3)
        test_GroupElement_interface(rand(F3, 2)...)
    end

end
