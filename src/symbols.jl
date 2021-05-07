change_pow(s::S, n::Integer) where S<:GSymbol = S(s.id, n)

function Base.iterate(s::GS, i=1) where GS<:GSymbol
    return i <= abs(s.pow) ? (change_pow(s, sign(s.pow)), i+1) : nothing
end
Base.size(s::GSymbol) = (abs(s.pow), )
Base.length(s::GSymbol) = first(size(s))

Base.eltype(s::GS) where GS<:GSymbol = GS

Base.isone(s::GSymbol) = iszero(s.pow)
Base.literal_pow(::typeof(^), s::Groups.GSymbol, ::Val{-1}) = inv(s)
Base.inv(s::GSymbol) = change_pow(s, -s.pow)
Base.hash(s::S, h::UInt) where S<:GSymbol = hash(s.id, hash(s.pow, hash(S, h)))

function Base.:(==)(s::GSymbol, t::GSymbol)
    isone(s) && isone(t) && return true
    s.pow == t.pow && s.id == t.id && return true
    return false
end

Base.convert(::Type{GS}, s::GSymbol) where GS<:GSymbol = GS(s.id, s.pow)
Base.convert(::Type{GS}, s::GS) where GS<:GSymbol = s
