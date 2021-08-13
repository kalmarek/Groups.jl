@testset "Aut(Σ₃.₀)" begin
    genus = 3

    π₁Σ = Groups.SurfaceGroup(genus, 0)

    Groups.PermRightAut(p::Perm) = Groups.PermRightAut(p.d)
    # Groups.PermLeftAut(p::Perm) = Groups.PermLeftAut(p.d)
    autπ₁Σ = let autπ₁Σ = AutomorphismGroup(π₁Σ)
        pauts = let p = perm"(1,3,5)(2,4,6)"
            [Groups.PermRightAut(p^i) for i in 0:2]
        end
        T = eltype(KnuthBendix.letters(alphabet(autπ₁Σ)))
        S = eltype(pauts)

        A = Alphabet(Union{T,S}[KnuthBendix.letters(alphabet(autπ₁Σ)); pauts])

        autG = AutomorphismGroup(
            π₁Σ,
            autπ₁Σ.gens,
            A,
            ntuple(i->inv(gens(π₁Σ, i)), 2Groups.genus(π₁Σ))
        )

        autG
    end

    Al = alphabet(autπ₁Σ)
    S = [gens(autπ₁Σ); inv.(gens(autπ₁Σ))]

    sautFn = let ltrs = KnuthBendix.letters(Al)
        parent(first(ltrs).autFn_word)
    end

    τ = Groups.rotation_element(sautFn)

    @testset "Twists" begin
        A = KnuthBendix.alphabet(sautFn)
        λ = Groups.ΡΛ(:λ, A, 2genus)
        ϱ = Groups.ΡΛ(:ϱ, A, 2genus)
        @test sautFn(Groups.Te_diagonal(λ, ϱ, 1)) ==
            conj(sautFn(Groups.Te_diagonal(λ, ϱ, 2)), τ)

        @test sautFn(Groups.Te_diagonal(λ, ϱ, 3)) == sautFn(Groups.Te(λ, ϱ, 3, 1))
    end

    z = let d = Groups.domain(τ)
        Groups.evaluate(τ^genus)
    end

    @test π₁Σ.(word.(z)) == Groups.domain(first(S))
    d = Groups.domain(first(S))
    p = perm"(1,3,5)(2,4,6)"
    @test Groups.evaluate!(deepcopy(d), τ) == d^inv(p)
    @test Groups.evaluate!(deepcopy(d), τ^2) == d^p

    E, sizes = Groups.wlmetric_ball(S, radius=3)
    @test sizes == [49, 1813, 62971]
    B2 = @view E[1:sizes[2]]

    σ = autπ₁Σ(Word([Al[Groups.PermRightAut(p)]]))

    @test conj(S[7], σ) == S[10]
    @test conj(S[7], σ^2) == S[11]
    @test conj(S[9], σ) == S[12]
    @test conj(S[9], σ^2) == S[8]

    @test conj(S[1], σ) == S[4]
    @test conj(S[1], σ^2) == S[5]
    @test conj(S[3], σ) == S[6]
    @test conj(S[3], σ^2) == S[2]

    B2ᶜ = [conj(b, σ) for b in B2]
    @test B2ᶜ != B2
    @test Set(B2ᶜ) == Set(B2)
end
