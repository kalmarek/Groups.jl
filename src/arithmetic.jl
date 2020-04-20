function Base.inv(W::T) where T<:GWord
    length(W) == 0 && return one(W)
    G = parent(W)
    w = T([inv(s) for s in Iterators.reverse(syllables(W))])
    return setparent!(w, G)
end

###############################################################################
#
#   Binary operators
#

function Base.push!(w::GWord{T}, s::T) where T <: GSymbol
    push!(syllables(w), s)
    return w
end

function Base.pushfirst!(w::GWord{T}, s::T) where T <: GSymbol
    pushfirst!(syllables(w), s)
    return w
end

function Base.append!(w::T, v::T) where T <: GWord
    append!(syllables(w), syllables(v))
    return w
end

function Base.prepend!(w::T, v::T) where T <: GWord
    prepend!(syllables(w), syllables(v))
    return w
end

Base.append!(w::T, v::T, others::Vararg{T,N}) where {N,T <: GWord} =
    append!(append!(w, v), others...)

function rmul!(out::T, x::T, y::T) where T<: GWord
    if out === x
        out = deepcopy(out)
        return freereduce!(append!(out, y))
    elseif out === y
        out = deepcopy(out)
        return freereduce!(prepend!(out, x))
    else
        slenx = syllablelength(x)
        sleny = syllablelength(y)
        resize!(syllables(out), slenx+sleny)
        syllables(out)[1:slenx] .= syllables(x)
        syllables(out)[slenx+1:slenx+sleny] .= syllables(y)
        return freereduce!(out)
    end
end

rmul!(out::T, v::T) where T<:GWord = freereduce!(append!(out, v))
lmul!(out::T, v::T) where T<:GWord = freereduce!(prepend!(out, v))

lmul!(out::T, x::T, y::T) where T <: GWord = rmul!(out, y, x)

AbstractAlgebra.mul!(out::T, x::T, y::T) where T <: GWord = rmul!(out, x, y)

(*)(W::GW, Z::GW) where GW <: GWord = rmul!(deepcopy(W), W, Z)
(*)(W::GWord, s::GSymbol) = freereduce!(push!(deepcopy(W), s))
(*)(s::GSymbol, W::GWord) = freereduce!(pushfirst!(deepcopy(W), s))

function power_by_squaring(W::GWord, p::Integer)
    if p < 0
        return power_by_squaring(inv(W), -p)
    elseif p == 0
        return one(W)
    elseif p == 1
        return W
    elseif p == 2
        return W*W
    end
    W = deepcopy(W)
    t = trailing_zeros(p) + 1
    p >>= t
    while (t -= 1) > 0
        append!(W, W)
    end
    Z = deepcopy(W)
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) >= 0
            append!(W, W)
        end
        append!(Z, W)
    end

    return freereduce!(Z)
end

(^)(x::GWord, n::Integer) = power_by_squaring(x,n)
