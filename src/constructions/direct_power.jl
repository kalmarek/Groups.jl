struct DirectPower{Gr,N,GEl<:GroupsCore.GroupElement} <: GroupsCore.Group
    group::Gr

    function DirectPower{N}(G::GroupsCore.Group) where {N}
        @assert N > 1
        return new{typeof(G),N,eltype(G)}(G)
    end
end

struct DirectPowerElement{GEl,N,Gr<:GroupsCore.Group} <: GroupsCore.GroupElement
    elts::NTuple{N,GEl}
    parent::DirectPower{Gr,N,GEl}
end

function DirectPowerElement(
    elts::AbstractVector{<:GroupsCore.GroupElement},
    G::DirectPower,
)
    return DirectPowerElement(ntuple(i -> elts[i], _nfold(G)), G)
end

_nfold(::DirectPower{Gr,N}) where {Gr,N} = N

function Base.one(G::DirectPower)
    return DirectPowerElement(ntuple(_ -> one(G.group), _nfold(G)), G)
end

function Base.eltype(::Type{<:DirectPower{Gr,N,GEl}}) where {Gr,N,GEl}
    return DirectPowerElement{GEl,N,Gr}
end

function Base.iterate(G::DirectPower)
    itr = Iterators.ProductIterator(ntuple(i -> G.group, _nfold(G)))
    res = iterate(itr)
    @assert res !== nothing
    elt = DirectPowerElement(first(res), G)
    return elt, (iterator = itr, state = last(res))
end

function Base.iterate(G::DirectPower, state)
    itr, st = state.iterator, state.state
    res = iterate(itr, st)
    res === nothing && return nothing
    elt = DirectPowerElement(first(res), G)
    return elt, (iterator = itr, state = last(res))
end

function Base.IteratorSize(::Type{<:DirectPower{Gr,N}}) where {Gr,N}
    Base.IteratorSize(Gr) isa Base.HasLength && return Base.HasShape{N}()
    Base.IteratorSize(Gr) isa Base.HasShape && return Base.HasShape{N}()
    return Base.IteratorSize(Gr)
end

Base.size(G::DirectPower) = ntuple(_ -> length(G.group), _nfold(G))

function GroupsCore.order(::Type{I}, G::DirectPower) where {I<:Integer}
    return convert(I, order(I, G.group)^_nfold(G))
end

GroupsCore.ngens(G::DirectPower) = _nfold(G) * ngens(G.group)

function GroupsCore.gens(G::DirectPower)
    N = _nfold(G)
    S = gens(G.group)
    tups = [ntuple(j -> (i == j ? s : one(s)), N) for i in 1:N for s in S]

    return [DirectPowerElement(elts, G) for elts in tups]
end

Base.isfinite(G::DirectPower) = isfinite(G.group)

GroupsCore.parent(g::DirectPowerElement) = g.parent

function Base.:(==)(g::DirectPowerElement, h::DirectPowerElement)
    return (parent(g) === parent(h) && g.elts == h.elts)
end

Base.hash(g::DirectPowerElement, h::UInt) = hash(g.elts, hash(parent(g), h))

Base.inv(g::DirectPowerElement) = DirectPowerElement(inv.(g.elts), parent(g))

function Base.:(*)(g::DirectPowerElement, h::DirectPowerElement)
    @assert parent(g) === parent(h)
    return DirectPowerElement(g.elts .* h.elts, parent(g))
end

# to make sure that parents are never copied i.e.
# g and deepcopy(g) share their parent
Base.deepcopy_internal(G::DirectPower, ::IdDict) = G

################## Implementing Group Interface Done!

function GroupsCore.gens(G::DirectPower, i::Integer)
    k = ngens(G.group)
    ci = CartesianIndices((k, _nfold(G)))
    @boundscheck checkbounds(ci, i)
    r, c = Tuple(ci[i])
    tup = ntuple(j -> j == c ? gens(G.group, r) : one(G.group), _nfold(G))
    return DirectPowerElement(tup, G)
end

# Overloading rand: the PRA of GroupsCore is known for not performing
# well on direct sums
function Random.Sampler(
    RNG::Type{<:Random.AbstractRNG},
    G::DirectPower,
    repetition::Random.Repetition = Val(Inf),
)
    return Random.SamplerTrivial(G)
end

function Base.rand(
    rng::Random.AbstractRNG,
    rs::Random.SamplerTrivial{<:DirectPower},
)
    G = rs[]
    return DirectPowerElement(rand(rng, G.group, _nfold(G)), G)
end

function GroupsCore.order(::Type{I}, g::DirectPowerElement) where {I<:Integer}
    return convert(I, reduce(lcm, (order(I, h) for h in g.elts); init = one(I)))
end

Base.isone(g::DirectPowerElement) = all(isone, g.elts)

function Base.show(io::IO, G::DirectPower)
    n = _nfold(G)
    nn = n == 1 ? "1-st" : n == 2 ? "2-nd" : n == 3 ? "3-rd" : "$n-th"
    return print(io, "Direct $(nn) power of ", G.group)
end

function Base.show(io::IO, g::DirectPowerElement)
    print(io, "( ")
    join(io, g.elts, ", ")
    return print(" )")
end

# convienience:
Base.@propagate_inbounds function Base.getindex(
    g::DirectPowerElement,
    i::Integer,
)
    return g.elts[i]
end
