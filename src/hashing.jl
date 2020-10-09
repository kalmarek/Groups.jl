###############################################################################
#
#   hashing, deepcopy and ==
#

function hash_internal(W::GWord)
    reduce!(W)
    h = hasparent(W) ? hash(parent(W)) : zero(UInt)
    return hash(syllables(W), hash(typeof(W), h))
end

function hash(W::GWord, h::UInt=UInt(0); kwargs...)
    if ismodified(W)
        savehash!(W, hash_internal(W; kwargs...))
        unsetmodified!(W)
    end
    return xor(savedhash(W), h)
end

# WARNING: Due to specialised (constant) hash function of GWords this one is actually necessary!
function Base.deepcopy_internal(W::T, dict::IdDict) where T<:GWord
    G = parent(W)
    g = T(deepcopy(syllables(W)))
    setparent!(g, G)
    return g
end

function (==)(W::T, Z::T) where T <: GWord
    hash(W) != hash(Z) && return false # distinguishes parent and parentless words
    if hasparent(W) && hasparent(Z)
        parent(W) != parent(Z) && return false
    end
    return syllables(W) == syllables(Z)
end
