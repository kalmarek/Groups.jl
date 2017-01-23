module Groups

import Base: length, ==, hash, show
import Base: one, inv, reduce, *, ^

export GSymbol, GWord

abstract GSymbol

function show(io::IO, s::GSymbol)
    if s.pow == 0
        print(io, "(id)")
    elseif s.pow == 1
        print(io, s.gen)
    else
        print(io, (s.gen)*"^$(s.pow)")
    end
end

length(s::GSymbol) = (s.pow == 0 ? 0 : 1)

IdSymbol(T::Type{GSymbol}) = throw(ArgumentError("Define IdSymbol(::Type{$T}) which is the identity element for Your type!"))

one{T<:GSymbol}(::Type{T}) = IdSymbol(T)
one(s::GSymbol) = one(typeof(s))

(*){T<:GSymbol}(s::T, t::T) = return GWord{T}([s])*t

change_pow(s::GSymbol, n::Int) = throw(ArgumentError("Define change_pow function for $(typeof(s))!"))


abstract Word

type GWord{T<:GSymbol} <: Word
    symbols::Vector{T}
    savedhash::UInt
    modified::Bool

    function GWord(symbols::Vector{T})
        return new(symbols, hash(symbols), true)
    end
end

GWord{T<:GSymbol}(s::T) = GWord{T}([s])

IDWord{T<:GSymbol}(::Type{T}) = GWord(one(T))
IDWord{T<:GSymbol}(W::GWord{T}) = IDWord(T)

function length(W::GWord)
    return sum([abs(s.pow) for s in W.symbols])
end

one{T}(::Type{GWord{T}}) = IDWord(T)
one{T}(w::GWord{T}) = one(GWord{T})

function inv{T}(W::GWord{T})
    if length(W) == 0
        return W
    else
        w = GWord{T}(reverse([inv(s) for s in W.symbols]))
        w.modified = true
        return w
    end
end

function join_free_symbols!(W::GWord)
    reduced = true
    for i in 1:length(W.symbols) - 1
        if W.symbols[i].gen == W.symbols[i+1].gen
            reduced = false
            p1 = W.symbols[i].pow
            p2 = W.symbols[i+1].pow
            W.symbols[i+1] = change_pow(W.symbols[i], p1 + p2)
            W.symbols[i] = one(W.symbols[i])
        end
    end
    return reduced
end

function freegroup_reduce!{T}(W::GWord{T})
    if length(W) < 2
        deleteat!(W.symbols, find(x -> x.pow == 0, W.symbols))
    else
        reduced = false
        while !reduced
            reduced = join_free_symbols!(W)
            deleteat!(W.symbols, find(x -> x.pow == 0, W.symbols))
        end
    end

    W.modified = false
    W.savedhash = hash(W.symbols,hash(typeof(W)))
    return W
end

freegroup_reduce(W::GWord) = freegroup_reduce!(deepcopy(W))

hash{T}(W::GWord{T}) = (W.modified && freegroup_reduce!(W); W.savedhash)

function (==){T}(W::GWord{T}, Z::GWord{T})
     W.modified && freegroup_reduce!(W) # reduce clears the flag and recalculate the hash
     Z.modified && freegroup_reduce!(Z)
     return W.savedhash == Z.savedhash && W.symbols == Z.symbols
end

(==){T}(W::GWord{T}, s::T) = W == GWord(s)
(==){T}(s::T, W::GWord{T}) = W == GWord(s)

function show(io::IO, W::GWord)
    if length(W) == 0
        print(io, "(id)")
    else
        join(io, [string(s) for s in W.symbols], "*")
    end
end

function r_multiply!(W::GWord, x; reduced::Bool=true)
    if length(x) > 0
        push!(W.symbols, x...)
    end
    if reduced
        freegroup_reduce!(W)
    end
    return W
end

function l_multiply!(W::GWord, x; reduced::Bool=true)
    if length(x) > 0
        unshift!(W.symbols, reverse(x)...)
    end
    if reduced
        freegroup_reduce!(W)
    end
    return W
end

r_multiply(W::GWord, x; reduced::Bool=true) =
    r_multiply!(deepcopy(W),x, reduced=reduced)
l_multiply(W::GWord, x; reduced::Bool=true) =
    l_multiply!(deepcopy(W),x, reduced=reduced)

(*){T}(W::GWord{T}, Z::GWord{T}) = r_multiply(W, Z.symbols)
(*)(W::GWord, s::GSymbol) = W*GWord(s)
(*)(s::GSymbol, W::GWord) = GWord(s)*W

function power_by_squaring{T}(x::GWord{T}, p::Integer)
    if p < 0
        return power_by_squaring(inv(x), -p)
    elseif p == 0
        return one(x)
    elseif p == 1
        return deepcopy(x)
    elseif p == 2
        return x*x
    end
    t = trailing_zeros(p) + 1
    p >>= t
    while (t -= 1) > 0
        x *= x
    end
    y = x
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) >= 0
            x *= x
        end
        y *= x
    end
    return freegroup_reduce!(y)
end

(^)(x::GWord, n::Integer) = power_by_squaring(x,n)
(^){T<:GSymbol}(x::T, n::Integer) = GWord(x)^n

include("free_groups.jl")

end # of module Groups
