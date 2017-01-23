import Base:convert

export FGSymbol, FGWord

immutable FGSymbol <: GSymbol
    gen::String
    pow::Int
end

(==)(s::FGSymbol, t::FGSymbol) = s.gen == t.gen && s.pow == t.pow
hash(s::FGSymbol, h::UInt) = hash(s.gen, hash(s.pow, hash(:FGSymbol, h)))

IdSymbol(::Type{FGSymbol}) = FGSymbol("(id)", 0)
FGSymbol(x::String) = FGSymbol(x,1)

inv(s::FGSymbol) = FGSymbol(s.gen, -s.pow)
convert(::Type{FGSymbol}, x::String) = FGSymbol(x)
change_pow(s::FGSymbol, n::Int) = (n==0 ? i=one(s) : FGSymbol(s.gen, n))

typealias FGWord GWord{FGSymbol}

FGWord(s::FGSymbol) = FGWord([s])
