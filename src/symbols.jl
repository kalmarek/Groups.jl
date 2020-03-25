change_pow(s::S, n::Integer) where S<:GSymbol = S(s.id, n)

function Base.iterate(s::GS, i=1) where GS<:GSymbol
    return i <= abs(s.pow) ? (GS(s.id, sign(s.pow)), i+1) : nothing
end
Base.length(s::GSymbol) = abs(s.pow)
Base.size(s::GSymbol) = (length(s), )
Base.eltype(s::GS) where GS<:GSymbol = GS

Base.isone(s::GSymbol) = iszero(s.pow)
Base.inv(s::GSymbol) = change_pow(s, -s.pow)
Base.hash(s::S, h::UInt) where S<:GSymbol = hash(s.id, hash(s.pow, hash(S, h)))

function (==)(s::GSymbol, t::GSymbol)
    isone(s) && isone(t) && return true
    s.pow == t.pow && s.id == t.id && return true
    return false
end
