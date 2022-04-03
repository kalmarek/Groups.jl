using PermutationGroups
using Groups.KnuthBendix

@testset "Wajnryb presentation for Σ₄₁" begin

    genus = 4

    Fn = FreeGroup(2genus)
    G = SpecialAutomorphismGroup(Fn)

    T = Groups.mcg_twists(G)

    # symplectic pairing in the free Group goes like this:
    # f1 ↔ f5
    # f2 ↔ f6
    # f3 ↔ f7
    # f4 ↔ f8

    T = let G = G
        (Tas, Tαs, Tes) = Groups.mcg_twists(G)
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
    a0 = (a1 * a2 * a3)^4 * b0^-1 # from the 3-chain relation
    X = a4 * a5 * a3 * a4 # auxillary, not present in the Primer
    b1 = X^-1 * a0 * X
    b2 = T[10]^-1 # Te₁₃

    As = T[[1, 5, 9, 6, 12, 7, 14, 8]] # the inverses of the elements a

    @testset "preserving relator" begin
        F = Groups.object(G)

        R = prod(commutator(gens(F,2i+1), gens(F,2i+2)) for i in 0:genus-1)

        for g in T
            @test g(R) == R
        end
    end

    @testset "commutation relations" begin
        for (i, ai) in enumerate(As) #the element ai here is actually the inverse of ai before. It does not matter for commutativity. Also, a0 is as defined before.
            for (j, aj) in enumerate(As)
                if abs(i - j) > 1
                    @test ai * aj == aj * ai
                elseif i ≠ j
                    @test ai * aj != aj * ai
                end
            end
            if i != 4
                @test a0 * ai == ai * a0
            end
        end
    end

    @testset "braid relations" begin
        for (i, ai) in enumerate(As) #the element ai here is actually the inverse of ai before. It does not matter for braid relations.
            for (j, aj) in enumerate(As)
                if abs(i - j) == 1
                    @test ai * aj * ai == aj * ai * aj
                end
            end
        end
        @test a0 * a4 * a0 == a4 * a0 * a4 # here, a0 and a4 are as before
    end

    @testset "3-chain relation" begin
        x = a4*a3*a2*a1*a1*a2*a3*a4 # auxillary; does not have a name in the Primer
        @test b0 == x*a0*x^-1
    end

    @testset "Lantern relation" begin

        @testset "b2 definition" begin
            @test b2 == (a2 * a3 * a1 * a2)^-1 * b1 * (a2 * a3 * a1 * a2)

            # some additional tests, checking what explicitly happens to the generators of the π₁ under b2
            d = Groups.domain(b2)
            img = evaluate(b2)
            z = img[3] * d[3]^-1

            @test img[1] == d[1]
            @test img[2] == d[2]
            @test img[3] == z * d[3]
            @test img[4] == z * d[4] * z^-1
            @test img[5] == z * d[5] * z^-1
            @test img[6] == z * d[6] * z^-1
            @test img[7] == d[7] * z^-1
            @test img[8] == d[8]
        end

        @testset "b2: commutation relations" begin
            @test b2 * a1 == a1 * b2
            @test b2 * a2 != a2 * b2
            @test b2 * a3 == a3 * b2
            @test b2 * a4 == a4 * b2
            @test b2 * a5 == a5 * b2
            @test b2 * a6 != a6 * b2
        end

        @testset "b2: braid relations" begin
            @test a2 * b2 * a2 == b2 * a2 * b2
            @test a6 * b2 * a6 == b2 * a6 * b2
        end

        @testset "lantern" begin
            u = (a6 * a5)^-1 * b1 * (a6 * a5)
            x = (a6 * a5 * a4 * a3 * a2 * u * a1^-1 * a2^-1 * a3^-1 * a4^-1) # yet another auxillary
            # x = (a4^-1*a3^-1*a2^-1*a1^-1*u*a2*a3*a4*a5*a6)

            b3 = x * a0 * x^-1
            b3im = evaluate(b3)
            b3cim = let g = b3
                f = Groups.compiled(g)
                f(Groups.domain(g))
            end
            @test b3im == b3cim
            @test a0 * b2 * b1 == a1 * a3 * a5 * b3

            @time evaluate(x)
            let g = b3
                f = Groups.compiled(g)
                f(Groups.domain(g))
                @time f(Groups.domain(g))
            end
        end
    end

    Base.conj(t::Groups.Transvection, p::Perm) =
        Groups.Transvection(t.id, t.i^p, t.j^p, t.inv)

    function Base.conj(elt::FPGroupElement, p::Perm)
        G = parent(elt)
        A = alphabet(elt)
        return G([A[conj(A[idx], p)] for idx in word(elt)])
    end

    @testset "Te₂₃ definition" begin
        Te₁₂, Te₂₃ = T[9], T[12]
        G = parent(Te₁₂)
        F₈ = Groups.object(G)
        (δ, d, γ, c, β, b, α, a) = Groups.gens(F₈)

        Groups.domain(Te₁₂)

        img_Te₂₃ = evaluate(Te₂₃)
        y = c * β^-1 * b^-1 * β
        @test img_Te₂₃ == (δ, d, y * γ, y * c * y^-1, β * y^-1, b, α, a)

        σ = perm"(7,5,3)(8,6,4)"
        Te₂₃_σ = conj(Te₁₂, σ)
        # @test word(Te₂₃_σ) == word(Te₂₃)

        @test evaluate(Te₂₃_σ) == evaluate(Te₂₃)
        @test Te₂₃ == Te₂₃_σ
    end

    @testset "Te₃₄ definition" begin
        Te₁₂, Te₃₄ = T[9], T[14]
        G = parent(Te₁₂)
        F₈ = Groups.object(G)
        (δ, d, γ, c, β, b, α, a) = Groups.gens(F₈)

        σ = perm"(7,3)(8,4)(5,1)(6,2)"
        Te₃₄_σ = conj(Te₁₂, σ)
        @test Te₃₄ == Te₃₄_σ
    end

    @testset "hyperelliptic τ is central" begin

        τ = Groups.rotation_element(G)
        τᵍ = τ^genus

        symplectic_gens = let genus = genus, G = G
            π₁Σ = Groups.SurfaceGroup(genus, 0)
            s_twists = Groups.symplectic_twists(π₁Σ)
            G.(word(t.autFn_word) for t in s_twists)
        end

        @test all(sg * τᵍ == τᵍ * sg for sg in symplectic_gens)
    end
end
