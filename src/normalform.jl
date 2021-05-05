"""
    normalform!(g::FPGroupElement)
Compute the normal form of `g`, possibly modifying `g` in-place.
"""
@inline function normalform!(g::FPGroupElement)
    isnormalform(g) && return g

    let w = one(word(g))
        w = normalform!(w, g)
        resize!(word(g), length(w))
        copyto!(word(g), w)
    end

    g = _update_savedhash!(g)
    @assert isnormalform(g)
    return g
end

"""
    normalform!(res::GEl, g::GEl) where GEl<:FPGroupElement
Compute the normal fom of `g`, storing it in `res`.
"""
function normalform!(res::GEl, g::GEl) where GEl<:FPGroupElement
    @boundscheck @assert parent(res) === parent(g)
    if isnormalform(g)
        copyto!(res, g)
        _update_savedhash!(res, g.savedhash)
    else
        resize!(word(res), 0)
        normalform!(word(res), g)
        res = _update_savedhash!(res)
    end
    return res
end

"""
    normalform!(res::AbstractWord, g::FPGroupElement)
Append the normal form of `g` to word `res`, modifying `res` in place.

Defaults to the rewriting in the free group.
"""
@inline function normalform!(res::AbstractWord, g::FPGroupElement)
    isone(res) && isnormalform(g) && return append!(res, word(g))
    if isnormalform(g) && inv(alphabet(g), last(out)) != first(word(g))
        return append!(res, word(g))
    end
    return KnuthBendix.rewrite_from_left!(res, word(g), rewriting(parent(g)))
end

"""
    free_rewrite!(v::AbstractWord, w::AbstractWord, A::Alphabet)
Append `w` to `v` applying free reductions as defined by the inverses of `A`.
"""
free_rewrite!(v::AbstractWord, w::AbstractWord, A::Alphabet) =
    KnuthBendix.rewrite_from_left!(v, w, A)

function KnuthBendix.rewrite_from_left!(
    v::AbstractWord,
    w::AbstractWord,
    A::Alphabet
    )
    while !isone(w)
        if isone(v)
            push!(v, popfirst!(w))
        else
            # the first check is for monoids only
            if KnuthBendix.hasinverse(last(v), A) && inv(A, last(v)) == first(w)
                pop!(v)
                popfirst!(w)
            else
                push!(v, popfirst!(w))
            end
        end
    end
    return v
end
