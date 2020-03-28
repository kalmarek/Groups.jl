function Base.inv(W::T) where T<:GWord
    length(W) == 0 && return W
    G = parent(W)
    w = T([inv(s) for s in Iterators.reverse(syllables(W))])
    return setparent!(w, G)
end

###############################################################################
#
#   Binary operators
#

function Base.append!(w::GWord{T}, v::AbstractVector{T}) where T
    append!(syllables(w), v)
    return w
end

function Base.prepend!(w::GWord{T}, v::AbstractVector{T}) where T
    prepend!(syllables(w), v)
    return w
end

Base.append!(w::T, v::T) where T <: GWord = append!(w, syllables(v))
Base.prepend!(w::T, v::T) where T <: GWord = prepend!(w, syllables(v))

for (mul, f) in ((:rmul!, :push!), (:lmul!, :pushfirst!))
    @eval begin
        function $mul(out::T, w::T, s::GSymbol) where T <:GWord
            resize!(syllables(out), syllablelength(w))
            syllables(out) .= syllables(w)
            $f(syllables(out), s)
            return freereduce!(out)
        end
    end
end

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

lmul!(out::T, x::T, y::T) where T <: GWord = rmul!(out, y, x)

function AbstractAlgebra.mul!(out::T, x::T, y::T) where T <: GWord
    return rmul!(out, x, y)
end

(*)(W::GW, Z::GW) where GW <: GWord = rmul!(deepcopy(W), W, Z)
(*)(W::GWord, s::GSymbol) = rmul!(deepcopy(W), W, s)
(*)(s::GSymbol, W::GWord) = lmul!(deepcopy(W), W, s)

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
