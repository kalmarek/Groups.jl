function test_homomorphism(hom)
    F = hom.source
    @test isone(hom(one(F)))
    @test all(inv(hom(g)) == hom(inv(g)) for g in gens(F))
    @test all(isone(hom(g)*hom(inv(g))) for g in gens(F))
    @test all(hom(g*h) == hom(g)*hom(h) for g in gens(F) for h in gens(F))
end

@testset "Homomorphisms" begin

    F₂ = FreeGroup(2)
    g,h = gens(F₂)

    ℤ² = FPGroup(F₂, [g*h => h*g])

    let hom = Groups.Homomorphism((i, G, H) -> Groups.word_type(H)([i]), F₂, ℤ²)

        @test hom(word(g)) == word(g)

        @test hom(word(g*h*inv(g))) == [1,3,2]

        @test hom(g*h*inv(g)) == hom(h)
        @test isone(hom(g*h*inv(g)*inv(h)))

        @test contains(sprint(print, hom), "Homomorphism")

        test_homomorphism(hom)
    end

    SAutF3 = SpecialAutomorphismGroup(FreeGroup(3))
    SL3Z = MatrixGroups.SpecialLinearGroup{3}(Int8)

    let hom = Groups.Homomorphism(
            Groups._abelianize,
            SAutF3,
            SL3Z,
        )

        A = alphabet(SAutF3)
        g = SAutF3([A[Groups.ϱ(1,2)]])
        h = SAutF3([A[Groups.λ(1,2)]])^-1

        @test !isone(g) && !isone(hom(g))
        @test !isone(h) && !isone(hom(h))
        @test !isone(g*h) && isone(hom(g*h))

        test_homomorphism(hom)
    end
end
