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

