import OrderedCollections: OrderedSet

mutable struct FPGroupIter{S,T,GEl}
    seen::S
    seen_iter_state::T
    current::GEl
    gen_idx::Int
    u_tmp::GEl
    v_tmp::GEl
end

function Base.iterate(G::AbstractFPGroup)
    seen = OrderedSet([one(G)])
    current, seen_state = iterate(seen)
    gr_iter = FPGroupIter(seen, seen_state, current, 1, one(G), one(G))
    return one(G), gr_iter
end

function _next_elt!(itr::FPGroupIter)
    res = iterate(itr.seen, itr.seen_iter_state)
    res === nothing && return true

    itr.current = first(res)
    itr.seen_iter_state = last(res)
    itr.gen_idx = 1

    return false
end

function Base.iterate(G::AbstractFPGroup, state)
    if state.gen_idx > length(alphabet(G))
        finished = _next_elt!(state)
        finished && return nothing
    end

    elt = let u = state.u_tmp, v = state.v_tmp, current = state.current
        copyto!(v, current)
        push!(word(v), state.gen_idx)
        _setnormalform!(v, false)
        _setvalidhash!(v, false)
        @assert !isnormalform(v)

        resize!(word(u), 0)
        normalform!(u, v)
    end

    state.gen_idx += 1

    if elt in state.seen
        return iterate(G, state)
    else
        @assert isnormalform(elt)
        push!(state.seen, deepcopy(elt))
        return elt, state
    end
end

# Groups.Core default:
# Base.IteratorSize(::Type{<:AbstractFPGroup}) = Base.SizeUnknown()
