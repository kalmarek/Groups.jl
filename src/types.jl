## "Abstract" definitions

"""
    AbstractFPGroup

An Abstract type representing finitely presented groups. Every instance `` must implement
 * `KnuthBendix.alphabet(G::MyFPGroup)`
 * `rewriting(G::MyFPGroup)` : return the rewriting object which must implement
 > `KnuthBendix.rewrite_from_left!(u, v, rewriting(G))`.
By default `alphabet(G)` is returned, which amounts to free rewriting in `G`.
 * `relations(G::MyFPGroup)` : return a set of defining relations.

AbstractFPGroup may also override `word_type(::Type{MyFPGroup}) = Word{UInt16}`,
which controls the word type used for group elements. If a group has more than `255` generators you need to define e.g.
> `word_type(::Type{MyFPGroup}) = Word{UInt16}`
"""
abstract type AbstractFPGroup <: GroupsCore.Group end

word_type(G::AbstractFPGroup) = word_type(typeof(G))
# the default:
word_type(::Type{<:AbstractFPGroup}) = Word{UInt8}

# the default (results in free rewriting)
rewriting(G::AbstractFPGroup) = alphabet(G)

Base.@propagate_inbounds function (G::AbstractFPGroup)(word::AbstractVector{<:Integer})
    @boundscheck @assert all(l -> 1 <= l <= length(KnuthBendix.alphabet(G)), word)
    return FPGroupElement(word_type(G)(word), G)
end

## Group Interface

Base.one(G::AbstractFPGroup) = FPGroupElement(one(word_type(G)), G)

Base.eltype(::Type{FPG}) where {FPG<:AbstractFPGroup} = FPGroupElement{FPG,word_type(FPG)}

include("iteration.jl")

GroupsCore.ngens(G::AbstractFPGroup) = length(G.gens)

function GroupsCore.gens(G::AbstractFPGroup, i::Integer)
    @boundscheck 1 <= i <= GroupsCore.ngens(G)
    l = alphabet(G)[G.gens[i]]
    return FPGroupElement(word_type(G)([l]), G)
end
GroupsCore.gens(G::AbstractFPGroup) = [gens(G, i) for i in 1:GroupsCore.ngens(G)]

# TODO: ProductReplacementAlgorithm
function Base.rand(rng::Random.AbstractRNG, rs::Random.SamplerTrivial{<:AbstractFPGroup})
    l = rand(10:100)
    G = rs[]
    nletters = length(alphabet(G))
    return FPGroupElement(word_type(G)(rand(1:nletters, l)), G)
end

Base.isfinite(::AbstractFPGroup) = (@warn "using generic isfinite(::AbstractFPGroup): the returned `false` might be wrong"; false)

## FPGroupElement

abstract type AbstractFPGroupElement{Gr} <: GroupElement end

mutable struct FPGroupElement{Gr<:AbstractFPGroup,W<:AbstractWord} <: AbstractFPGroupElement{Gr}
    word::W
    savedhash::UInt
    parent::Gr

    FPGroupElement(word::W, G::AbstractFPGroup, hash::UInt=UInt(0)) where {W<:AbstractWord} =
        new{typeof(G),W}(word, hash, G)

    FPGroupElement{Gr, W}(word::AbstractWord, G::Gr) where {Gr, W} =
        new{Gr, W}(word, UInt(0), G)
end

word(f::AbstractFPGroupElement) = f.word

#convenience
KnuthBendix.alphabet(g::AbstractFPGroupElement) = alphabet(parent(g))

function Base.show(io::IO, f::AbstractFPGroupElement)
    f = normalform!(f)
    KnuthBendix.print_repr(io, word(f), alphabet(f))
end

## GroupElement Interface for FPGroupElement

Base.parent(f::AbstractFPGroupElement) = f.parent

function Base.:(==)(g::AbstractFPGroupElement, h::AbstractFPGroupElement)
    @boundscheck @assert parent(g) === parent(h)
    normalform!(g)
    normalform!(h)
    hash(g) != hash(h) && return false
    return word(g) == word(h)
end

function Base.deepcopy_internal(g::FPGroupElement, stackdict::IdDict)
    return FPGroupElement(copy(word(g)), parent(g), g.savedhash)
end

function Base.inv(g::GEl) where GEl <: AbstractFPGroupElement
    G = parent(g)
    return GEl(inv(alphabet(G), word(g)), G)
end

function Base.:(*)(g::GEl, h::GEl) where GEl<:AbstractFPGroupElement
    @boundscheck @assert parent(g) === parent(h)
    return GEl(word(g) * word(h), parent(g))
end

GroupsCore.isfiniteorder(g::AbstractFPGroupElement) = isone(g) ? true : (@warn "using generic isfiniteorder(::AbstractFPGroupElement): the returned `false` might be wrong"; false)

# additional methods:
Base.isone(g::AbstractFPGroupElement) = (normalform!(g); isempty(word(g)))

## Free Groups

struct FreeGroup{T} <: AbstractFPGroup
    gens::Vector{T}
    alphabet::KnuthBendix.Alphabet{T}

    function FreeGroup(gens, A::KnuthBendix.Alphabet) where {W}
        @assert length(gens) == length(unique(gens))
        @assert all(l -> l in KnuthBendix.letters(A), gens)
        return new{eltype(gens)}(gens, A)
    end
end

function FreeGroup(A::Alphabet)
    @boundscheck @assert all(KnuthBendix.hasinverse(l, A) for l in KnuthBendix.letters(A))
    ltrs = KnuthBendix.letters(A)
    gens = Vector{eltype(ltrs)}()
    invs = Vector{eltype(ltrs)}()
    for l in ltrs
        l ∈ invs && continue
        push!(gens, l)
        push!(invs, inv(A, l))
    end

    return FreeGroup(gens, A)
end

function FreeGroup(n::Integer)
    symbols = Symbol[]
    inverses = Int[]
    sizehint!(symbols, 2n)
    sizehint!(inverses, 2n)
    for i in 1:n
        push!(symbols, Symbol(:f, i), Symbol(:F, i))
        push!(inverses, 2i, 2i-1)
    end
    return FreeGroup(symbols[1:2:2n], Alphabet(symbols, inverses))
end

Base.show(io::IO, F::FreeGroup) = print(io, "free group on $(ngens(F)) generators")

# mandatory methods:
KnuthBendix.alphabet(F::FreeGroup) = F.alphabet
relations(F::FreeGroup) = Pair{eltype(F)}[]

# GroupsCore interface:
# these are mathematically correct
Base.isfinite(::FreeGroup) = false

GroupsCore.isfiniteorder(g::AbstractFPGroupElement{<:FreeGroup}) = isone(g) ? true : false

## FP Groups

struct FPGroup{T,R,S} <: AbstractFPGroup
    gens::Vector{T}
    relations::Vector{Pair{S,S}}
    rws::R
end

KnuthBendix.alphabet(G::FPGroup) = alphabet(rewriting(G))
rewriting(G::FPGroup) = G.rws

relations(G::FPGroup) = G.relations

function FPGroup(
    G::AbstractFPGroup,
    rels::AbstractVector{<:Pair{GEl,GEl}};
    ordering = KnuthBendix.LenLex,
    kwargs...,
) where {GEl<:FPGroupElement}

    O = ordering(alphabet(G))
    for (lhs, rhs) in rels
        @assert parent(lhs) === parent(rhs) === G
    end
    word_rels = [word(lhs) => word(rhs) for (lhs, rhs) in [relations(G); rels]]
    rws = KnuthBendix.RewritingSystem(word_rels, O)

    KnuthBendix.knuthbendix!(rws; kwargs...)

    return FPGroup(G.gens, rels, rws)
end

function Base.show(io::IO, G::FPGroup)
    print(io, "⟨")
    join(io, gens(G), ", ")
    print(io, " | ")
    join(io, relations(G), ", ")
    print(io, "⟩")
end

## GSymbol aka letter of alphabet

abstract type GSymbol end
Base.literal_pow(::typeof(^), t::GSymbol, ::Val{-1}) = inv(t)

function subscriptify(n::Integer)
    subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
    return join([Char(subscript_0 + i) for i in reverse(digits(n))], "")
end
