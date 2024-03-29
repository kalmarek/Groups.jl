## "Abstract" definitions

"""
    AbstractFPGroup

An Abstract type representing finitely presented groups. Every instance must implement
 * `KnuthBendix.alphabet(G::MyFPGroup)`
 * `rewriting(G::MyFPGroup)` : return the rewriting object which must implement
 > `KnuthBendix.rewrite!(u, v, rewriting(G))`.
 E.g. for `G::FreeGroup` `alphabet(G)` is returned, which amounts to free rewriting.
 * `ordering(G::MyFPGroup)[ = KnuthBendix.ordering(rewriting(G))]` : return the
 (implicit) ordering for the alphabet of `G`.
 * `relations(G::MyFPGroup)` : return a set of defining relations.

AbstractFPGroup may also override `word_type(::Type{MyFPGroup}) = Word{UInt8}`,
which controls the word type used for group elements.
If a group has more than `255` generators you need to define e.g.
> `word_type(::Type{MyFPGroup}) = Word{UInt16}`
"""
abstract type AbstractFPGroup <: GroupsCore.Group end

word_type(G::AbstractFPGroup) = word_type(typeof(G))
# the default:
word_type(::Type{<:AbstractFPGroup}) = Word{UInt8}

"""
    rewriting(G::AbstractFPGroup)
Return a "rewriting object" for elements of `G`.

The rewriting object must must implement
    KnuthBendix.rewrite!(u::AbstractWord, v::AbstractWord, rewriting(G))

For example if `G` is a `FreeGroup` then `alphabet(G)` is returned which results
in free rewriting. For `FPGroup` a rewriting system is returned which may
(or may not) rewrite word `v` to its normal form (depending on e.g. its confluence).
"""
function rewriting end

KnuthBendix.ordering(G::AbstractFPGroup) = ordering(rewriting(G))
KnuthBendix.alphabet(G::AbstractFPGroup) = alphabet(ordering(G))

Base.@propagate_inbounds function (G::AbstractFPGroup)(
    word::AbstractVector{<:Integer},
)
    @boundscheck @assert all(l -> 1 <= l <= length(alphabet(G)), word)
    return FPGroupElement(word_type(G)(word), G)
end

## Group Interface

Base.one(G::AbstractFPGroup) = FPGroupElement(one(word_type(G)), G)

function Base.eltype(::Type{FPG}) where {FPG<:AbstractFPGroup}
    return FPGroupElement{FPG,word_type(FPG)}
end

include("iteration.jl")

GroupsCore.ngens(G::AbstractFPGroup) = length(G.gens)

function GroupsCore.gens(G::AbstractFPGroup, i::Integer)
    @boundscheck 1 <= i <= GroupsCore.ngens(G)
    l = alphabet(G)[G.gens[i]]
    return FPGroupElement(word_type(G)([l]), G)
end
function GroupsCore.gens(G::AbstractFPGroup)
    return [gens(G, i) for i in 1:GroupsCore.ngens(G)]
end

function Base.isfinite(::AbstractFPGroup)
    return (
        @warn "using generic isfinite(::AbstractFPGroup): the returned `false` might be wrong"; false
    )
end

## FPGroupElement

abstract type AbstractFPGroupElement{Gr} <: GroupElement end

Base.copy(g::AbstractFPGroupElement) = one(g) * g
word(f::AbstractFPGroupElement) = f.word

mutable struct FPGroupElement{Gr<:AbstractFPGroup,W<:AbstractWord} <:
               AbstractFPGroupElement{Gr}
    word::W
    savedhash::UInt
    parent::Gr

    function FPGroupElement(
        word::W,
        G::AbstractFPGroup,
        hash::UInt = UInt(0),
    ) where {W<:AbstractWord}
        return new{typeof(G),W}(word, hash, G)
    end

    function FPGroupElement{Gr,W}(word::AbstractWord, G::Gr) where {Gr,W}
        return new{Gr,W}(word, UInt(0), G)
    end
end

function Base.copy(f::FPGroupElement)
    return FPGroupElement(copy(word(f)), parent(f), f.savedhash)
end

#convenience
KnuthBendix.alphabet(g::AbstractFPGroupElement) = alphabet(parent(g))

function Base.show(io::IO, f::AbstractFPGroupElement)
    f = normalform!(f)
    return KnuthBendix.print_repr(io, word(f), alphabet(f))
end

## GroupElement Interface for FPGroupElement

Base.parent(f::AbstractFPGroupElement) = f.parent

function Base.:(==)(g::AbstractFPGroupElement, h::AbstractFPGroupElement)
    @boundscheck @assert parent(g) === parent(h)
    normalform!(g)
    normalform!(h)
    # I. compare hashes of the normalform
    # II. compare some data associated to FPGroupElement,
    # e.g. word, image of the domain etc.
    hash(g) != hash(h) && return false
    equality_data(g) == equality_data(h) && return true # compares

    # if this failed it is still possible that the words together can be
    # rewritten even further, so we
    # 1. rewrite word(g⁻¹·h) w.r.t. rewriting(parent(g))
    # 2. check if the result is empty
    G = parent(g)

    g⁻¹h = append!(inv(word(g), alphabet(G)), word(h))
    # similar + empty preserve the storage size
    # saves some re-allocations if res does not represent id
    res = similar(word(g))
    resize!(res, 0)
    res = KnuthBendix.rewrite!(res, g⁻¹h, rewriting(G))
    return isone(res)
end

function Base.deepcopy_internal(g::FPGroupElement, stackdict::IdDict)
    haskey(stackdict, g) && return stackdict[g]
    cw = Base.deepcopy_internal(word(g), stackdict)
    h = FPGroupElement(cw, parent(g), g.savedhash)
    stackdict[g] = h
    return h
end

function Base.inv(g::GEl) where {GEl<:AbstractFPGroupElement}
    G = parent(g)
    return GEl(inv(word(g), alphabet(G)), G)
end

function Base.:(*)(g::GEl, h::GEl) where {GEl<:AbstractFPGroupElement}
    @boundscheck @assert parent(g) === parent(h)
    A = alphabet(parent(g))
    k = 0
    while k + 1 ≤ min(length(word(g)), length(word(h)))
        if inv(word(g)[end-k], A) == word(h)[k+1]
            k += 1
        else
            break
        end
    end
    w = @view(word(g)[1:end-k]) * @view(word(h)[k+1:end])
    res = GEl(w, parent(g))
    return res
end

function GroupsCore.isfiniteorder(g::AbstractFPGroupElement)
    return isone(g) ? true :
           (
        @warn "using generic isfiniteorder(::AbstractFPGroupElement): the returned `false` might be wrong"; false
    )
end

# additional methods:
Base.isone(g::AbstractFPGroupElement) = (normalform!(g); isempty(word(g)))

## Free Groups

struct FreeGroup{T,O} <: AbstractFPGroup
    gens::Vector{T}
    ordering::O

    function FreeGroup(gens, ordering::KnuthBendix.WordOrdering)
        @assert length(gens) == length(unique(gens))
        @assert all(l -> l in alphabet(ordering), gens)
        return new{eltype(gens),typeof(ordering)}(gens, ordering)
    end
end

function FreeGroup(n::Integer)
    symbols =
        collect(Iterators.flatten((Symbol(:f, i), Symbol(:F, i)) for i in 1:n))
    inverses = collect(Iterators.flatten((2i, 2i - 1) for i in 1:n))
    return FreeGroup(Alphabet(symbols, inverses))
end

FreeGroup(A::Alphabet) = FreeGroup(KnuthBendix.LenLex(A))

function __group_gens(A::Alphabet)
    @boundscheck @assert all(KnuthBendix.hasinverse(l, A) for l in A)
    gens = Vector{eltype(A)}()
    invs = Vector{eltype(A)}()
    for l in A
        l ∈ invs && continue
        push!(gens, l)
        push!(invs, inv(l, A))
    end
    return gens
end

function FreeGroup(O::KnuthBendix.WordOrdering)
    grp_gens = __group_gens(alphabet(O))
    return FreeGroup(grp_gens, O)
end

function Base.show(io::IO, F::FreeGroup)
    return print(io, "free group on $(ngens(F)) generators")
end

# mandatory methods:
KnuthBendix.ordering(F::FreeGroup) = F.ordering
rewriting(F::FreeGroup) = alphabet(F) # alphabet(F) = alphabet(ordering(F))
relations(F::FreeGroup) = Pair{eltype(F),eltype(F)}[]

# GroupsCore interface:
# these are mathematically correct
Base.isfinite(::FreeGroup) = false

function GroupsCore.isfiniteorder(g::AbstractFPGroupElement{<:FreeGroup})
    return isone(g) ? true : false
end

## FP Groups

struct FPGroup{T,RW,S} <: AbstractFPGroup
    gens::Vector{T}
    relations::Vector{Pair{S,S}}
    rw::RW
end

relations(G::FPGroup) = G.relations
rewriting(G::FPGroup) = G.rw

function FPGroup(
    G::AbstractFPGroup,
    rels::AbstractVector{<:Pair{GEl,GEl}};
    ordering = KnuthBendix.ordering(G),
    kwargs...,
) where {GEl<:FPGroupElement}
    for (lhs, rhs) in rels
        @assert parent(lhs) === parent(rhs) === G
    end
    word_rels = [word(lhs) => word(rhs) for (lhs, rhs) in [relations(G); rels]]
    rws = KnuthBendix.RewritingSystem(word_rels, ordering)

    rws = KnuthBendix.knuthbendix(rws, KnuthBendix.Settings(; kwargs...))

    return FPGroup(G.gens, rels, KnuthBendix.IndexAutomaton(rws))
end

function Base.show(io::IO, ::MIME"text/plain", G::FPGroup)
    println(
        io,
        "Finitely presented group generated by $(ngens(G)) element",
        ngens(G) > 1 ? 's' : "",
        ": ",
    )
    join(io, gens(G), ", ", ", and ")
    println(
        io,
        "\n    subject to relation",
        length(relations(G)) > 1 ? 's' : "",
    )
    return Base.print_array(io, relations(G))
end

function Base.show(io::IO, G::FPGroup)
    print(io, "⟨")
    Base.print_array(io, permutedims(gens(G)))
    println(io, " | ")
    print(io, "\t ")
    Base.print_array(io, permutedims(relations(G)))
    return print(io, " ⟩")
end

function Base.show(io::IO, ::Type{<:FPGroup{T}}) where {T}
    return print(io, FPGroup, "{$T, …}")
end

## GSymbol aka letter of alphabet

abstract type GSymbol end
Base.literal_pow(::typeof(^), t::GSymbol, ::Val{-1}) = inv(t)

function subscriptify(n::Integer)
    subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
    return join([Char(subscript_0 + i) for i in reverse(digits(n))], "")
end
