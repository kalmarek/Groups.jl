import PermutationGroups as PG

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
struct WreathProduct{DP<:DirectPower,PGr<:PG.AbstractPermutationGroup} <:
       GroupsCore.Group
    N::DP
    P::PGr

    function WreathProduct(G::Group, P::PG.AbstractPermutationGroup)
        N = DirectPower{PG.AP.degree(P)}(G)
        return new{typeof(N),typeof(P)}(N, P)
    end
end

struct WreathProductElement{
    DPEl<:DirectPowerElement,
    PEl<:PG.AP.AbstractPermutation,
    Wr<:WreathProduct,
} <: GroupsCore.GroupElement
    n::DPEl
    p::PEl
    parent::Wr

    function WreathProductElement(
        n::DirectPowerElement,
        p::PG.AP.AbstractPermutation,
        W::WreathProduct,
    )
        return new{typeof(n),typeof(p),typeof(W)}(n, p, W)
    end
end

Base.one(W::WreathProduct) = WreathProductElement(one(W.N), one(W.P), W)

function Base.eltype(::Type{<:WreathProduct{DP,PGr}}) where {DP,PGr}
    return WreathProductElement{eltype(DP),eltype(PGr),WreathProduct{DP,PGr}}
end

function Base.iterate(G::WreathProduct)
    itr = Iterators.product(G.N, G.P)
    res = iterate(itr)
    @assert res !== nothing
    ab, st = res
    (a, b) = ab
    elt = WreathProductElement(a, b, G)
    return elt, (itr, st)
end

function Base.iterate(G::WreathProduct, state)
    itr, st = state
    res = iterate(itr, st)
    res === nothing && return nothing
    (a::eltype(G.N), b::eltype(G.P)), st = res
    elt = WreathProductElement(a, b, G)
    return elt, (itr, st)
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

function GroupsCore.order(::Type{I}, G::WreathProduct) where {I<:Integer}
    return convert(I, order(I, G.N) * order(I, G.P))
end

function GroupsCore.gens(G::WreathProduct)
    N_gens = [WreathProductElement(n, one(G.P), G) for n in gens(G.N)]
    P_gens = [WreathProductElement(one(G.N), p, G) for p in gens(G.P)]
    return [N_gens; P_gens]
end

Base.isfinite(G::WreathProduct) = isfinite(G.N) && isfinite(G.P)

GroupsCore.parent(g::WreathProductElement) = g.parent

function Base.:(==)(g::WreathProductElement, h::WreathProductElement)
    return parent(g) === parent(h) && g.n == h.n && g.p == h.p
end

function Base.hash(g::WreathProductElement, h::UInt)
    return hash(g.n, hash(g.p, hash(g.parent, h)))
end

function _act(p::PG.AP.AbstractPermutation, n::DirectPowerElement)
    return DirectPowerElement(
        ntuple(i -> n.elts[i^p], length(n.elts)),
        parent(n),
    )
end

function Base.inv(g::WreathProductElement)
    pinv = inv(g.p)
    return WreathProductElement(_act(pinv, inv(g.n)), pinv, parent(g))
end

function Base.:(*)(g::WreathProductElement, h::WreathProductElement)
    @assert parent(g) === parent(h)
    return WreathProductElement(g.n * _act(g.p, h.n), g.p * h.p, parent(g))
end

# to make sure that parents are never copied i.e.
# g and deepcopy(g) share their parent
Base.deepcopy_internal(G::WreathProduct, ::IdDict) = G

################## Implementing Group Interface Done!

# Overloading rand: the PRA of GroupsCore is known for not performing
# well on direct sums
function Random.Sampler(
    RNG::Type{<:Random.AbstractRNG},
    G::WreathProduct,
    repetition::Random.Repetition = Val(Inf),
)
    return Random.SamplerTrivial(G)
end

function Base.rand(
    rng::Random.AbstractRNG,
    rs::Random.SamplerTrivial{<:WreathProduct},
)
    G = rs[]
    return WreathProductElement(rand(rng, G.N), rand(rng, G.P), G)
end

Base.isone(g::WreathProductElement) = isone(g.n) && isone(g.p)

function Base.show(io::IO, G::WreathProduct)
    return print(io, "Wreath product of ", G.N.group, " by ", G.P)
end

function Base.show(io::IO, g::WreathProductElement)
    return print(io, "( ", g.n, "≀", g.p, " )")
end
