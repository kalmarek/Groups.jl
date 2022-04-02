struct DirectProduct{Gt,Ht,GEl,HEl} <: GroupsCore.Group
    first::Gt
    last::Ht

    function DirectProduct(G::GroupsCore.Group, H::GroupsCore.Group)
        return new{typeof(G),typeof(H),eltype(G),eltype(H)}(G, H)
    end
end

struct DirectProductElement{GEl,HEl,Gt,Ht} <: GroupsCore.GroupElement
    elts::Tuple{GEl,HEl}
    parent::DirectProduct{Gt,Ht,GEl,HEl}
end

DirectProductElement(g, h, G::DirectProduct) = DirectProduct((g, h), G)

Base.one(G::DirectProduct) =
    DirectProductElement((one(G.first), one(G.last)), G)

Base.eltype(::Type{<:DirectProduct{Gt,Ht,GEl,HEl}}) where {Gt,Ht,GEl,HEl} =
    DirectProductElement{GEl,HEl,Gt,Ht}

function Base.iterate(G::DirectProduct)
    itr = Iterators.product(G.first, G.last)
    res = iterate(itr)
    @assert res !== nothing
    elt = DirectProductElement(first(res), G)
    return elt, (iterator = itr, state = last(res))
end

function Base.iterate(G::DirectProduct, state)
    itr, st = state.iterator, state.state
    res = iterate(itr, st)
    res === nothing && return nothing
    elt = DirectProductElement(first(res), G)
    return elt, (iterator = itr, state = last(res))
end

function Base.IteratorSize(::Type{<:DirectProduct{Gt,Ht}}) where {Gt,Ht}
    Gi = Base.IteratorSize(Gt)
    Hi = Base.IteratorSize(Ht)
    if Gi isa Base.IsInfinite || Hi isa Base.IsInfinite
        return Base.IsInfinite()
    elseif Gi isa Base.SizeUnknown || Hi isa Base.SizeUnknown
        return Base.SizeUnknown()
    else
        return Base.HasShape{2}()
    end
end

Base.size(G::DirectProduct) = (length(G.first), length(G.last))

GroupsCore.order(::Type{I}, G::DirectProduct) where {I<:Integer} =
    convert(I, order(I, G.first) * order(I, G.last))

GroupsCore.ngens(G::DirectProduct) = ngens(G.first) + ngens(G.last)

function GroupsCore.gens(G::DirectProduct)
    gens_first = [DirectProductElement((g, one(G.last)), G) for g in gens(G.first)]

    gens_last = [DirectProductElement((one(G.first), g), G) for g in gens(G.last)]

    return [gens_first; gens_last]
end

Base.isfinite(G::DirectProduct) = isfinite(G.first) && isfinite(G.last)

function Base.rand(
    rng::Random.AbstractRNG,
    rs::Random.SamplerTrivial{<:DirectProduct},
)
    G = rs[]
    return DirectProductElement((rand(rng, G.first), rand(rng, G.last)), G)
end

GroupsCore.parent(g::DirectProductElement) = g.parent

Base.:(==)(g::DirectProductElement, h::DirectProductElement) =
    (parent(g) === parent(h) && g.elts == h.elts)

Base.hash(g::DirectProductElement, h::UInt) = hash(g.elts, hash(parent(g), h))

Base.deepcopy_internal(g::DirectProductElement, stackdict::IdDict) =
    DirectProductElement(Base.deepcopy_internal(g.elts, stackdict), parent(g))

Base.inv(g::DirectProductElement) =
    DirectProductElement(inv.(g.elts), parent(g))

function Base.:(*)(g::DirectProductElement, h::DirectProductElement)
    @assert parent(g) === parent(h)
    return DirectProductElement(g.elts .* h.elts, parent(g))
end

GroupsCore.order(::Type{I}, g::DirectProductElement) where {I<:Integer} =
    convert(I, lcm(order(I, first(g.elts)), order(I, last(g.elts))))

Base.isone(g::DirectProductElement) = all(isone, g.elts)

Base.show(io::IO, G::DirectProduct) =
    print(io, "Direct product of $(G.first) and $(G.last)")
Base.show(io::IO, g::DirectProductElement) =
    print(io, "( $(join(g.elts, ",")) )")
