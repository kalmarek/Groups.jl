function _abelianize(
    i::Integer,
    source::AutomorphismGroup{<:FreeGroup},
    target::MatrixGroups.SpecialLinearGroup{N, T}) where {N, T}
    n = ngens(object(source))
    @assert n == N
    aut = alphabet(source)[i]
    if aut isa Transvection
        # we change (i,j) to (j, i) to be consistent with the action:
        # Automorphisms act on the right which corresponds to action on
        # the columns in the matrix case
        eij = MatrixGroups.ElementaryMatrix{N}(
                aut.j,
                aut.i,
                ifelse(aut.inv, -one(T), one(T))
            )
        k = alphabet(target)[eij]
        return word_type(target)([k])
    else
        throw("unexpected automorphism symbol: $(typeof(aut))")
    end
end

function _abelianize(
    i::Integer,
    source::AutomorphismGroup{<:Groups.SurfaceGroup},
    target::MatrixGroups.SpecialLinearGroup{N, T}) where {N, T}
    n = ngens(Groups.object(source))
    @assert n == N
    g = alphabet(source)[i].autFn_word
    result = one(target)
    for l in word(g)
        append!(word(result), _abelianize(l, parent(g), target))
    end

    return word(result)
end

function Groups._abelianize(
    i::Integer,
    source::AutomorphismGroup{<:Groups.SurfaceGroup},
    target::MatrixGroups.SymplecticGroup{N, T}
) where {N, T}
    @assert iseven(N)
    As = alphabet(source)
    At = alphabet(target)

    SlN = let genus = Groups.genus(Groups.object(source))
        @assert 2genus == N
        MatrixGroups.SpecialLinearGroup{2genus}(T)
    end

    ab = Groups.Homomorphism(Groups._abelianize, source, SlN, check=false)

    matrix_spn_map = let S = gens(target)
        Dict(MatrixGroups.matrix_repr(g)=> word(g) for g in union(S, inv.(S)))
    end

    # renumeration:
    # (f1, f2, f3, f4, f5, f6) = (a₁, a₂, a₃, b₁, b₂, b₃) →
    # → (b₃, a₃, b₂, a₂, b₁, a₁)
    # hence p = [6, 4, 2, 5, 3, 1]
    p = [reverse(2:2:N); reverse(1:2:N)]

    g = source([i])
    Mg = MatrixGroups.matrix_repr(ab(g))[p,p]

    return matrix_spn_map[Mg]
end
