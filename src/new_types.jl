using GroupsCore
# using Groups
# import Groups.AbstractFPGroup
import KnuthBendix
import KnuthBendix: AbstractWord, Alphabet, Word, RewritingSystem
import KnuthBendix: alphabet
using Random

## "Abstract" definitions

"""
    AbstractFPGroup

An Abstract type representing finitely presented groups. Every instance `` must implement
 * `KnuthBendix.alphabet(G::MyFPGroup)`
 * `rewriting(G::MyFPGroup)` : return the rewriting object which must implement
 > `KnuthBendix.rewrite_from_left!(u, v, rewriting(G))`.
By default `alphabet(G)` is returned, which amounts to free rewriting in `G`.

AbstractFPGroup may also override `word_type(::Type{MyFPGroup}) = Word{UInt16}`,
which controls the word type used for group elements. If if your group has less
than `255` generators you may define
> `word_type(::Type{MyFPGroup}) = Word{UInt8}`
"""
abstract type AbstractFPGroup <: GroupsCore.Group end

word_type(G::AbstractFPGroup) = word_type(typeof(G))
# the default:
word_type(::Type{<:AbstractFPGroup}) = Word{UInt16}

rewriting(G::AbstractFPGroup) = alphabet(G)

function (G::AbstractFPGroup)(word::AbstractVector{<:Integer})
    @boundscheck @assert all(l -> 1<= l <=length(KnuthBendix.alphabet(G)), word)
    return FPGroupElement(word_type(G)(word), G)
end

## Group Interface

Base.one(G::AbstractFPGroup) = FPGroupElement(one(word_type(G)), G)

Base.eltype(::Type{FPG}) where {FPG<:AbstractFPGroup} =
    FPGroupElement{FPG, word_type(FPG)}

struct FPGroupIter{GEl}
    elts::Vector{GEl}
    seen::Set{GEl}
    u::GEl
    v::GEl
end

FPGroupIter(G::AbstractFPGroup) =
    FPGroupIter([one(G)], Set([one(G)]), one(G), one(G))

Base.iterate(G::AbstractFPGroup) = one(G), (FPGroupIter(G), 1, 1)
@inline function Base.iterate(G::AbstractFPGroup, state)
    iter, elt_idx, gen_idx = state

    if gen_idx > length(alphabet(G))
        elt_idx == length(iter.elts) && return nothing
        gen_idx = 1
        elt_idx += 1
    end


    res = let (u, v) = (iter.u, iter.v), elt = iter.elts[elt_idx]
        copyto!(v, elt) # this invalidates normalform of v
        @assert !isnormalform(v)
        push!(word(v), gen_idx)
        resize!(word(u), 0)

        normalform!(u, v)
    end

    if res in iter.seen
        return iterate(G, (iter, elt_idx, gen_idx+1))
    else
        w = deepcopy(res)
        @assert isnormalform(w)
        push!(iter.elts, w)
        push!(iter.seen, w)
        state = (iter, elt_idx, gen_idx+1)
        return w, state
    end
end

# the default:
# Base.IteratorSize(::Type{<:AbstractFPGroup}) = Base.SizeUnknown()

GroupsCore.ngens(G::AbstractFPGroup) = length(G.gens)

function GroupsCore.gens(G::AbstractFPGroup, i::Integer)
    @boundscheck 1<=i<=GroupsCore.ngens(G)
    l = alphabet(G)[G.gens[i]]
    return FPGroupElement(word_type(G)([l]), G)
end
GroupsCore.gens(G::AbstractFPGroup) = [gens(G, i) for i in 1:GroupsCore.ngens(G)]

# TODO: ProductReplacementAlgorithm
function Base.rand(
    rng::Random.AbstractRNG,
    rs::Random.SamplerTrivial{<:AbstractFPGroup},
    )
    l = rand(10:100)
    G = rs[]
    nletters = length(alphabet(G))
    return FPGroupElement(word_type(G)(rand(1:nletters, l)), G)
end

## FPGroupElement

mutable struct FPGroupElement{G<:AbstractFPGroup, W<:AbstractWord} <: GroupElement
    word::W
    savedhash::UInt
    parent::G

    FPGroupElement(word::W, G::AbstractFPGroup) where W<:AbstractWord =
        new{typeof(G), W}(word, UInt(0), G)

    FPGroupElement(word::W, hash::UInt, G::AbstractFPGroup) where W<:AbstractWord =
        new{typeof(G), W}(word, hash, G)
end

word(f::FPGroupElement) = f.word

#convenience
KnuthBendix.alphabet(g::FPGroupElement) = alphabet(parent(g))

function Base.show(io::IO, f::FPGroupElement)
    f = normalform!(f)
    print(io, KnuthBendix.string_repr(word(f), alphabet(f)))
end

## Hashing

# We store hash of a word in field `savedhash` to use it as cheap replacement
# for equality checking. In the last bit of `savehash` we store the information
# whether the word is already reduced (in normal form) or not. Hence
isnormalform(g::FPGroupElement) = Bool(g.savedhash & 1)

# To update hash use this internal method, possibly only after computing the
# normal form of `g`:
_update_savedhash!(g::FPGroupElement, h = hash(word(g))) =
    (g.savedhash = h | 1; return g)

# If `g` has been mutated (e.g. `word(g)` has been modified, `_invalidate_hash`
# must be called.
_invalidate_hash!(g::FPGroupElement) = g.savedhash = g.savedhash & 0

# Accessor for the saved hash
@inline _savedhash(f::FPGroupElement) = f.savedhash

function Base.hash(g::FPGroupElement, h::UInt)
    normalform!(g)
    # compute the best approximation of the normal form
    return hash(parent(g), _savedhash(g) ⊻ h)
end

function Base.copyto!(res::FPGroupElement, g::FPGroupElement)
    resize!(word(res), length(word(g)))
    copyto!(word(res), word(g))
    _invalidate_hash!(res)
    return res
end

## GroupElement Interface for FPGroupElement

Base.parent(f::FPGroupElement) = f.parent
GroupsCore.parent_type(::Type{<:FPGroupElement{G}}) where G = G

function Base.:(==)(g::FPGroupElement, h::FPGroupElement)
    @boundscheck @assert parent(g) === parent(h)
    hash(g) != hash(h) && return false
    # hash reduces to normal form behind the scenes
    return word(g) == word(h)
end

function Base.deepcopy_internal(g::FPGroupElement, stackdict::IdDict)
    return FPGroupElement(copy(word(g)), _savedhash(g), parent(g))
end

Base.inv(g::FPGroupElement) =
    (G = parent(g); FPGroupElement(inv(alphabet(G), word(g)), G))

function Base.:(*)(g::FPGroupElement, h::FPGroupElement)
    @boundscheck @assert parent(g) === parent(h)
    return FPGroupElement(word(g)*word(h), parent(g))
end

GroupsCore.isfiniteorder(g::FPGroupElement) = isone(g) ? true : throw("Not Implemented")

# additional methods:
Base.isone(g::FPGroupElement) = (normalform!(g); isempty(word(g)))

## Free Groups

struct FreeGroup{T} <: AbstractFPGroup
    gens::Vector{T}
    alphabet::KnuthBendix.Alphabet{T}

    function FreeGroup(gens, A::KnuthBendix.Alphabet) where W
        @assert length(gens) == length(unique(gens))
        @assert all(l->l in KnuthBendix.letters(A), gens)
        return new{eltype(gens)}(gens, A)
    end
end

function FreeGroup(A::Alphabet)
    @boundscheck @assert all(KnuthBendix.hasinverse(l, A)
        for l in KnuthBendix.letters(A))
    return FreeGroup(KnuthBendix.letters(A), A)
end

Base.show(io::IO, F::FreeGroup) = print(io, "free group on $(ngens(F)) generators")

# mandatory methods:
KnuthBendix.alphabet(F::FreeGroup) = F.alphabet
