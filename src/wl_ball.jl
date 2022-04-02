"""
    wlmetric_ball(S::AbstractVector{<:GroupElem}
        [, center=one(first(S)); radius=2, op=*])
Compute metric ball as a list of elements of non-decreasing length, given the
word-length metric on the group generated by `S`. The ball is centered at `center`
(by default: the identity element). `radius` and `op` keywords specify the
radius and multiplication operation to be used.
"""
function wlmetric_ball_serial(S::AbstractVector{T}, center::T=one(first(S)); radius = 2, op = *) where {T}
    @assert radius >= 1
    old = union!([center], [center*s for s in S])
    return _wlmetric_ball(S, old, radius, op, collect, unique!)
end

function wlmetric_ball_thr(S::AbstractVector{T}, center::T=one(first(S)); radius = 2, op = *) where {T}
    @assert radius >= 1
    old = union!([center], [center*s for s in S])
    return _wlmetric_ball(S, old, radius, op, Folds.collect, Folds.unique)
end

function _wlmetric_ball(S, old, radius, op, collect, unique)
    sizes = [1, length(old)]
    for r in 2:radius
        old = let old = old, S = S,
            new = collect(
                (g = op(o, s); hash(g); g)
                for o in @view(old[sizes[end-1]:end]) for s in S
            )
            append!(old, new)
            unique(old)
        end
        push!(sizes, length(old))
    end
    return old, sizes[2:end]
end

function wlmetric_ball(
    S::AbstractVector{T},
    center::T = one(first(S));
    radius = 2,
    op = *,
    threading = true,
) where {T}
    threading && return wlmetric_ball_thr(S, center, radius = radius, op = op)
    return wlmetric_ball_serial(S, center, radius = radius, op = op)
end
