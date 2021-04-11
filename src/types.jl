abstract type AbstractFPGroup <: GroupsCore.Group end

function Base.one(G::Gr) where Gr <: AbstractFPGroup
    El = eltype(G)
    id = El(eltype(El)[])
    id.parent = G
    return id
end

"""
    GSymbol
Represents a syllable. Abstract type which all group symbols of
`AbstractFPGroups` should subtype. Each concrete subtype should implement fields:
 * `id` which is the `Symbol` representation/identification of a symbol
 * `pow` which is the (multiplicative) exponent of a symbol.
"""
abstract type GSymbol end

abstract type GWord{T<:GSymbol} <: GroupsCore.GroupElement end

"""
    W::GroupWord{T} <: GWord{T<:GSymbol} <:GroupElem
Basic representation of element of a finitely presented group.
* `syllables(W)` return particular group syllables which multiplied constitute `W`
group as a word in generators.
* `parent(W)` return the parent group.

As the reduction (inside the parent group) of word to normal form may be time
consuming, we provide a shortcut that is useful in practice:
`savehash!(W, h)` and `ismodified(W)` functions.
When computing `hash(W)`, a reduction to normal form is performed and a
persistent hash is stored inside `W`, setting `ismodified(W)` flag to `false`.
This hash can be accessed by `savedhash(W)`.
Future comparisons of `W` try not to perform reduction and use the stored hash as shortcut. Only when hashes collide reduction is performed. Whenever word `W` is
changed, `ismodified(W)` returns `false` and stored hash is invalidated.
"""

mutable struct GroupWord{T} <: GWord{T}
    symbols::Vector{T}
    modified::Bool
    savedhash::UInt
    parent::Group

    function GroupWord{T}(symbols::AbstractVector{<:GSymbol}) where T
       return new{T}(symbols, true, zero(UInt))
    end
    GroupWord(v::AbstractVector{T}) where T<:GSymbol = GroupWord{T}(v)
    GroupWord{T}(s::GSymbol) where T<:GSymbol = GroupWord{T}(T[s])
    GroupWord(s::T) where T<:GSymbol = GroupWord{T}(s)
end
