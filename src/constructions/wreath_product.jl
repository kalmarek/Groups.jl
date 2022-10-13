import PermutationGroups: AbstractPermutationGroup, AbstractPerm, degree, SymmetricGroup

"""
    WreathProduct(G::Group, P::AbstractPermutationGroup) <: Group
Return the wreath product of a group `G` by permutation group `P`, usually
written as `G ≀ P`.

As set `G ≀ P` is the same as `Gᵈ × P` and the group can be understood as a
semi-direct product of `P` acting on `d`-fold cartesian product of `G` by
permuting coordinates. To be more precise, the multiplication inside wreath
product is defined as
>  `(n, σ) * (m, τ) = (n*(m^σ), σ*τ)`
where `m^σ` denotes the action (from the right) of the permutation `σ` on
`d`-tuples of elements from `G`.
"""
struct WreathProduct{DP<:DirectPower,PGr<:AbstractPermutationGroup} <:
       GroupsCore.Group
    N::DP
    P::PGr

    function WreathProduct(G::Group, P::AbstractPermutationGroup)
        N = DirectPower{degree(P)}(G)
        return new{typeof(N),typeof(P)}(N, P)
    end
end

struct WreathProductElement{
    DPEl<:DirectPowerElement,
    PEl<:AbstractPerm,
    Wr<:WreathProduct,
} <: GroupsCore.GroupElement
    n::DPEl
    p::PEl
    parent::Wr

    function WreathProductElement(
        n::DirectPowerElement,
        p::AbstractPerm,
        W::WreathProduct,
    )
        new{typeof(n),typeof(p),typeof(W)}(n, p, W)
    end
end

Base.one(W::WreathProduct) = WreathProductElement(one(W.N), one(W.P), W)

Base.eltype(::Type{<:WreathProduct{DP,PGr}}) where {DP,PGr} =
    WreathProductElement{eltype(DP),eltype(PGr),WreathProduct{DP,PGr}}

function Base.iterate(G::WreathProduct)
    itr = Iterators.product(G.N, G.P)
    res = iterate(itr)
    @assert res !== nothing
    elt = WreathProductElement(first(res)..., G)
    return elt, (iterator=itr, state=last(res))
end

function Base.iterate(G::WreathProduct, state)
    itr, st = state.iterator, state.state
    res = iterate(itr, st)
    res === nothing && return nothing
    elt = WreathProductElement(first(res)..., G)
    return elt, (iterator=itr, state=last(res))
end

function Base.IteratorSize(::Type{<:WreathProduct{DP,PGr}}) where {DP,PGr}
    dpI = Base.IteratorSize(DP)
    pgI = Base.IteratorSize(PGr)

    if dpI isa Base.IsInfinite || pgI isa Base.IsInfinite
        return Base.IsInfinite()
    elseif dpI isa Base.SizeUnknown || pgI isa Base.SizeUnknown
        return Base.SizeUnknown()
    else
        return Base.HasShape{2}()
    end
end

Base.size(G::WreathProduct) = (length(G.N), length(G.P))

GroupsCore.order(::Type{I}, G::WreathProduct) where {I<:Integer} =
    convert(I, order(I, G.N) * order(I, G.P))

function GroupsCore.gens(G::WreathProduct)
    N_gens = [WreathProductElement(n, one(G.P), G) for n in gens(G.N)]
    P_gens = [WreathProductElement(one(G.N), p, G) for p in gens(G.P)]
    return [N_gens; P_gens]
end

Base.isfinite(G::WreathProduct) = isfinite(G.N) && isfinite(G.P)

function Base.rand(
    rng::Random.AbstractRNG,
    rs::Random.SamplerTrivial{<:WreathProduct},
)

    G = rs[]
    return WreathProductElement(rand(rng, G.N), rand(rng, G.P), G)
end

GroupsCore.parent(g::WreathProductElement) = g.parent

Base.:(==)(g::WreathProductElement, h::WreathProductElement) =
    parent(g) === parent(h) && g.n == h.n && g.p == h.p

Base.hash(g::WreathProductElement, h::UInt) =
    hash(g.n, hash(g.p, hash(g.parent, h)))

function Base.deepcopy_internal(g::WreathProductElement, stackdict::IdDict)
    return WreathProductElement(
        Base.deepcopy_internal(g.n, stackdict),
        Base.deepcopy_internal(g.p, stackdict),
        parent(g),
    )
end

_act(p::AbstractPerm, n::DirectPowerElement) =
    DirectPowerElement(n.elts^p, parent(n))

function Base.inv(g::WreathProductElement)
    pinv = inv(g.p)
    return WreathProductElement(_act(pinv, inv(g.n)), pinv, parent(g))
end

function Base.:(*)(g::WreathProductElement, h::WreathProductElement)
    @assert parent(g) === parent(h)
    return WreathProductElement(g.n * _act(g.p, h.n), g.p * h.p, parent(g))
end

Base.isone(g::WreathProductElement) = isone(g.n) && isone(g.p)

Base.show(io::IO, G::WreathProduct) =
    print(io, "Wreath product of $(G.N.group) by $(G.P)")
Base.show(io::IO, g::WreathProductElement) = print(io, "( $(g.n)≀$(g.p) )")

Base.copy(g::WreathProductElement) = WreathProductElement(g.n, g.p, parent(g))
