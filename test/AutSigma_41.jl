using PermutationGroups

@testset "Wajnryb presentation for Σ₄₁" begin

    genus = 4

    G = New.SpecialAutomorphismGroup(New.FreeGroup(2genus))

    T = let G = G; (Tas, Tαs, Tes) = New.mcg_twists(genus)
        Ta = G.(Tas)
        Tα = G.(Tαs)
        Tes = G.(Tes)

        [Ta; Tα; Tes]
    end

    a1 = T[1]^-1 # Ta₁
    a2 = T[5]^-1 # Tα₁
    a3 = T[9]^-1 # Te₁₂
    a4 = T[6]^-1 # Tα₂
    a5 = T[12]^-1 # Te₂₃
    a6 = T[7]^-1 # Tα₃
    a7 = T[14]^-1 # Te₃₄
    a8 = T[8]^-1 # Tα₄

    b0 = T[2]^-1 # Ta₂
    a0 = (a1*a2*a3)^4*b0^-1 # from the 3-chain relation
    X = a4*a5*a3*a4 # auxillary, not present in the Primer
    b1 = X^-1*a0*X
    b2 = T[10]^-1 # Te₁₃

    As = T[[1,5,9,6,12,7,14,8]] # the inverses of the elements a

    @testset "commutation relations" begin
        for (i, ai) in enumerate(As) #the element ai here is actually the inverse of ai before. It does not matter for commutativity. Also, a0 is as defined before.
            for (j, aj) in enumerate(As)
                if abs(i-j) > 1
                    @test ai*aj == aj*ai
                elseif i ≠ j
                    @test ai*aj != aj*ai
                end
            end
            if i != 4
                @test a0*ai == ai*a0
            end
        end
    end

    @testset "braid relations" begin
        for (i, ai) in enumerate(As) #the element ai here is actually the inverse of ai before. It does not matter for braid relations.
            for (j, aj) in enumerate(As)
                if abs(i-j) == 1
                    @test ai*aj*ai == aj*ai*aj
                end
            end
        end
        @test a0*a4*a0 == a4*a0*a4 # here, a0 and a4 are as before
    end

    @testset "Lantern relation" begin

        @testset "b2 definition" begin
                @test b2 == (a2*a3*a1*a2)^-1*b1*(a2*a3*a1*a2)

            # some additional tests, checking what explicitly happens to the generators of the π₁ under b2
            d = New.domain(b2)
            im = New.evaluate(b2)
            z = im[7]*d[7]^-1

            @test im[1] == d[1]
            @test im[2] == z*d[2]*z^-1
            @test im[3] == z*d[3]*z^-1
            @test im[4] == d[4]

            @test im[5] == d[5]*z^-1
            @test im[6] == z*d[6]*z^-1
            @test im[7] == z*d[7]
            @test im[8] == d[8]

        end

        @testset "b2: commutation relations" begin
            @test b2*a1 == a1*b2
            @test b2*a2 != a2*b2
            @test b2*a3 == a3*b2
            @test b2*a4 == a4*b2
            @test b2*a5 == a5*b2
            @test b2*a6 != a6*b2
        end

        @testset "b2: braid relations" begin
            @test a2*b2*a2 == b2*a2*b2
            @test a6*b2*a6 == b2*a6*b2
        end

        @testset "lantern" begin
            u = (a6*a5)^-1*b1*(a6*a5)
            x = (a6*a5*a4*a3*a2*u*a1^-1*a2^-1*a3^-1*a4^-1) # yet another auxillary
            # x = (a4^-1*a3^-1*a2^-1*a1^-1*u*a2*a3*a4*a5*a6)
            @time New.evaluate(x)
            b3 = x*a0*x^-1
            @time New.evaluate(b3)
            @test a0*b2*b1 == a1*a3*a5*b3
        end
    end


    @testset "Te₁₂ definition" begin
        G = parent(first(T))
        F₈ = New.object(G)
        (a, b, c, d, α, β, γ, δ) = gens(F₈)

        A = KnuthBendix.alphabet(G)

        λ = [i == j ? one(G) : G([A[New.λ(i,j)]]) for i in 1:8, j in 1:8]
        ϱ = [i == j ? one(G) : G([A[New.ϱ(i,j)]]) for i in 1:8, j in 1:8]

        g = one(G)
        # @show g
        # @show g(Groups.domain(G))

        # β ↦ α*β
        g *= λ[6,5]
        @test New.evaluate(g)[6] == α*β

        # α ↦ a*α*b^-1
        g *= λ[5,1]*inv(ϱ[5,2])
        @test New.evaluate(g)[5] == a*α*b^-1

        # β ↦ b*α^-1*a^-1*α*β
        g *= inv(λ[6,5])
        @test New.evaluate(g)[6] == b*α^-1*a^-1*α*β

        # b ↦ α
        g *= λ[2,5]*inv(λ[2,1]);
        @test New.evaluate(g)[2] == α

        # b ↦ b*α^-1*a^-1*α
        g *= inv(λ[2,5]);
        @test New.evaluate(g)[2] == b*α^-1*a^-1*α

        # b ↦ b*α^-1*a^-1*α*b*α^-1
        g *= inv(ϱ[2,5])*ϱ[2,1];
        @test New.evaluate(g)[2] == b*α^-1*a^-1*α*b*α^-1

        # b ↦ b*α^-1*a^-1*α*b*α^-1*a*α*b^-1
        g *= ϱ[2,5];
        @test New.evaluate(g)[2] == b*α^-1*a^-1*α*b*α^-1*a*α*b^-1

        x = b*α^-1*a^-1*α
        @test New.evaluate(g) == # (a, b,        c, d, α,      β,   γ, δ)
                                   (a, x*b*x^-1, c, d, α*x^-1, x*β, γ, δ)
        @test g == T[9]
    end

    Base.conj(t::New.Transvection, p::Perm) =
        New.Transvection(t.id, t.i^p, t.j^p, t.inv)

    function Base.conj(elt::New.FPGroupElement, p::Perm)
        G = parent(elt)
        A = New.alphabet(elt)
        return G([A[conj(A[idx], p)] for idx in New.word(elt)])
    end

    @testset "Te₂₃ definition" begin
        Te₁₂, Te₂₃ = T[9], T[12]
        G = parent(Te₁₂)
        F₈ = New.object(G)
        (a, b, c, d, α, β, γ, δ) = gens(F₈)

        img_Te₂₃ = New.evaluate(Te₂₃)
        y = c*β^-1*b^-1*β
        @test img_Te₂₃ == (a, b, y*c*y^-1, d, α, β*y^-1, y*γ, δ)

        σ = perm"(1,2,3)(5,6,7)(8)"
        Te₂₃_σ = conj(Te₁₂, σ)
        # @test New.word(Te₂₃_σ) == New.word(Te₂₃)

        @test New.evaluate(Te₂₃_σ) == New.evaluate(Te₂₃)
        @test Te₂₃ == Te₂₃_σ
    end

    @testset "Te₃₄ definition" begin
        Te₁₂, Te₃₄ = T[9], T[14]
        G = parent(Te₁₂)
        F₈ = New.object(G)
        (a, b, c, d, α, β, γ, δ) = Groups.gens(F₈)

        σ = perm"(1,3)(2,4)(5,7)(6,8)"
        Te₃₄_σ = conj(Te₁₂, σ)
        @test Te₃₄ == Te₃₄_σ
    end
end
