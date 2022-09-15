# Code by Marek Kaluba

# using OrderedCollections

mutable struct Subgroup{Gr, GrE} <: GroupsCore.Group
    supergroup::Gr
    gens::Vector{GrE}
    order::Int
    elts::OrderedSet{GrE}

    function Subgroup(G::Group, gens::AbstractVector{<:GroupElement};
        full_set=false)
        Gr = typeof(G)
        GrE = eltype(gens)
        H = new{Gr, GrE}(G, gens)
        if full_set
            @assert one(G) ∈ gens throw(ArgumentError("Full set of elements defining subgroup does not contain the identity."))
            H.order = length(gens)
            H.elts = OrderedSet(gens)
        end
        return H
    end
end

Base.parent(sg::Subgroup) = sg.supergroup

struct SubgroupElement{GrE, Gr<:Subgroup} <: GroupsCore.GroupElement
    elt::GrE
    parent::Gr
end

Base.one(sg::Subgroup) = SubgroupElement(one(parent(sg)), sg)
Base.eltype(::Type{<:Subgroup{Gr, GrE}}) where {Gr, GrE} =
    SubgroupElement{GrE, Subgroup{Gr, GrE}}

# may be false for finite subgroups of infinite groups
function Base.IteratorSize(::Type{<:Subgroup{Gr}}) where Gr
    PIS = Base.IteratorSize(Gr)
    if PIS in (Base.IsInfinite(), Base.SizeUnknown())
        return PIS
    else
        return Base.HasLength()
    end
end

###### iteration and order
mutable struct GroupIter{S, T, GEl}
    seen::S
    seen_iter_state::T
    current::GEl
    gen_idx::Int
    tmp::GEl
end

function _iterate(G::Subgroup)
    seen = OrderedSet([one(G)])
    current, seen_state = iterate(seen)
    gr_iter = GroupIter(seen, seen_state, current, 1, one(G))
    return one(G), gr_iter
end

function _next_elt!(itr::GroupIter)
    # iterate over itr.seen, updating itr.seen_iter_state
    # we're advancing to the next element, so
    # * set itr.current and
    # * reset gen_idx to 1

    res = iterate(itr.seen, itr.seen_iter_state)
    res == nothing && return true
    # returned true terminates the outer iteration

    itr.current = first(res)
    itr.seen_iter_state = last(res)
    itr.gen_idx = 1

    return false
end

function _iterate(G::Subgroup, state)

    # we generate new elements of G as elt*s^±1, where s is a generator;
    # gen_idx ranges from 1 to 2ngens(G), even idx denote inverses

    if state.gen_idx > 2ngens(G) # we've finished the generators, so advance current to the next element
        finished = _next_elt!(state)
        finished && return nothing
    end
    elt = let s = gens(G, (state.gen_idx + 1) ÷ 2)
        # multiply current by s (or its inverse) on the right
        state.tmp = GroupsCore.mul!(state.tmp, state.current, (isodd(state.gen_idx) ? s : inv(s)))
    end

    state.gen_idx += 1

    if elt in state.seen
        iterate(G, state)
    else
        push!(state.seen, deepcopy(elt))
        return deepcopy(elt), state
    end
end

function GroupsCore.order(::Type{T}, sg::Subgroup) where T
    if !isdefined(sg, :order)
    #morally: T(length(collect(sg)))
        g, gr_iter = _iterate(sg)

        while true
            k = _iterate(sg, gr_iter)
            isnothing(k) && break
        end
        sg.order = length(gr_iter.seen)
        sg.elts = gt_iter.seen
    end
    return convert(T, sg.order)
end

function Base.iterate(G::Subgroup)
    if !isdefined(G, :order)
        k = order(Int, G)
    end
    @assert isdefined(G, :elts)
    g, state = iterate(G.elts)
    return SubgroupElement(g, G), state
end

function Base.iterate(G::Subgroup, state)
    k = iterate(G.elts, state)
    isnothing(k) && return nothing
    g, state = k
    return SubgroupElement(g, G), state
end
###### iteration and order

GroupsCore.gens(sg::Subgroup) = [SubgroupElement(g, sg) for g in sg.gens]
GroupsCore.gens(sg::Subgroup, i::Integer) = SubgroupElement(sg.gens[i], sg)

function Base.rand(
    rng::Random.AbstractRNG,
    rs::Random.SamplerTrivial{<:Subgroup},
)
    sg = rs[]
    S = sg.gens
    @warn "FIXME: distribution of random element is not guaranteed to be uniform"
    return SubgroupElement(prod(rand(S, 10*length(S))), sg)
end

Base.parent(g::SubgroupElement) = g.parent

Base.:(==)(g::SubgroupElement, h::SubgroupElement) =
    (parent(g) === parent(h)) && g.elt == h.elt

Base.hash(g::SubgroupElement, h::UInt) = hash(g.elt, hash(parent(g), h))

function Base.deepcopy_internal(g::SubgroupElement, stackdict::IdDict)
    haskey(stackdict, g) && return stackdict[g]
    return SubgroupElement(Base.deepcopy_internal(g.elt, stackdict), parent(g))
end

# TODO:
Base.copy(g::SubgroupElement) = deepcopy(g)

Base.inv(g::SubgroupElement) = SubgroupElement(inv(g.elt), parent(g))

function Base.:(*)(g::SubgroupElement, h::SubgroupElement)
    @assert parent(g) === parent(h)
    return SubgroupElement(g.elt * h.elt, parent(g))
end

Base.show(io::IO, sg::Subgroup) =
    print(io, "subgroup of $(parent(sg)) defined by $(ngens(sg)) generators.")
Base.show(io::IO, g::SubgroupElement) = show(io, g.elt)
