"""
    normalform!(g::FPGroupElement)
Compute the normal form of `g`, possibly modifying `g` in-place.
"""
@inline function normalform!(g::AbstractFPGroupElement; force = false)
    if !force
        isnormalform(g) && return g
    end

    let w = one(word(g))
        w = normalform!(w, g)
        resize!(word(g), length(w))
        copyto!(word(g), w)
    end

    _setnormalform!(g, true)
    _setvalidhash!(g, false)
    @assert isnormalform(g)
    return g
end

"""
    normalform!(res::GEl, g::GEl) where GEl<:FPGroupElement
Compute the normal fom of `g`, storing it in `res`.
"""
function normalform!(res::GEl, g::GEl) where {GEl<:AbstractFPGroupElement}
    @boundscheck @assert parent(res) === parent(g)
    if isnormalform(g)
        copyto!(res, g)
    else
        resize!(word(res), 0)
        normalform!(word(res), g)
        _setnormalform!(res, true)
        _setvalidhash!(res, false)
    end
    return res
end

"""
    normalform!(res::AbstractWord, g::FPGroupElement)
Write the normal form of `g` to word `res`, modifying `res` in place.

The particular implementation of the normal form depends on `parent(g)`.
"""
@inline function normalform!(res::AbstractWord, g::AbstractFPGroupElement)
    isone(res) && isnormalform(g) && return append!(res, word(g))
    return KnuthBendix.rewrite!(res, word(g), rewriting(parent(g)))
end
