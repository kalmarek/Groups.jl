mutable struct FPIterState{GEl,T}
    elts::OrderedSet{GEl}
    u::GEl
    v::GEl
    elts_iter_state::T
end

Base.in(x, itr::FPIterState) = x in itr.elts
Base.push!(itr::FPIterState, x) = push!(itr.elts, x)
function FPIterState(G::AbstractFPGroup)
    S = OrderedSet([one(G)])
    elt, state = iterate(S)
    return elt, FPIterState(S, one(G), one(G), state)
end

function Base.iterate(G::AbstractFPGroup)
    elt, state = FPIterState(G)
    return one(G), (state, elt, 1)
end

function Base.iterate(G::AbstractFPGroup, state)
    iter, elt, gen_idx = state

    if gen_idx > length(alphabet(G))
        res = iterate(iter.elts, iter.elts_iter_state)
        res === nothing && return nothing
        gen_idx = 1
        elt = first(res)
        iter.elts_iter_state = last(res)
    end

    res = let (u, v) = (iter.u, iter.v), elt = elt
        copyto!(v, elt) # this invalidates normalform of v
        push!(word(v), gen_idx)
        resize!(word(u), 0)

        _setnormalform!(v, false)
        _setvalidhash!(v, false)

        @assert !isnormalform(v)
        normalform!(u, v)
    end

    if res in iter
        return iterate(G, (iter, elt, gen_idx + 1))
    else
        w = deepcopy(res)
        @assert isnormalform(w)
        push!(iter, w)
        return w, (iter, elt, gen_idx + 1)
    end
end

# Groups.Core default:
# Base.IteratorSize(::Type{<:AbstractFPGroup}) = Base.SizeUnknown()
