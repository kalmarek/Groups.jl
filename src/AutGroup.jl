###############################################################################
#
#   AutSymbol/ AutGroup / Automorphism
#
###############################################################################

struct RTransvect
    i::Int
    j::Int
end

struct LTransvect
    i::Int
    j::Int
end

struct FlipAut
    i::Int
end

struct PermAut
    perm::Nemo.Generic.perm{Int8}
end

struct Identity end

struct AutSymbol <: GSymbol
   str::String
   pow::Int
   typ::Union{LTransvect, RTransvect, PermAut, FlipAut, Identity}
end

mutable struct AutGroup{N} <: AbstractFPGroup
   objectGroup::FreeGroup
   gens::Vector{AutSymbol}
end

mutable struct Automorphism{N} <: GWord{AutSymbol}
    symbols::Vector{AutSymbol}
    modified::Bool
    savedhash::UInt
    parent::AutGroup{N}

    Automorphism{N}(f::Vector{AutSymbol}) where N = new(f, true)

end

export Automorphism, AutGroup

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

elem_type(::AutGroup{N}) where N = Automorphism{N}

parent_type(::Automorphism{N}) where N = AutGroup{N}

###############################################################################
#
#   AutSymbol defining functions
#
###############################################################################

function (ϱ::RTransvect)(v, pow=1::Int)
    @inbounds Groups.r_multiply!(v[ϱ.i], (v[ϱ.j]^pow).symbols, reduced=false)
    return v
end

function (λ::LTransvect)(v, pow=1::Int)
    @inbounds Groups.l_multiply!(v[λ.i], (v[λ.j]^pow).symbols, reduced=false)
    return v
end

function (σ::PermAut)(v, pow=1::Int)
   w = deepcopy(v)
   s = (σ.perm^pow).d
   @inbounds for k in eachindex(v)
       v[k].symbols = w[s[k]].symbols
   end
   return v
end

function (ɛ::FlipAut)(v, pow=1::Int)
   @inbounds if isodd(pow)
       v[ɛ.i].symbols = inv(v[ɛ.i]).symbols
   end
   return v
end

(::Identity)(v, pow=1::Int) = v

# taken from ValidatedNumerics, under under the MIT "Expat" License:
# https://github.com/JuliaIntervals/ValidatedNumerics.jl/blob/master/LICENSE.md
function subscriptify(n::Integer)
    subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
    return join([Char(subscript_0 + i) for i in reverse(digits(n))])
end

function id_autsymbol()
   return AutSymbol("(id)", 0, Identity())
end

function rmul_autsymbol(i, j; pow::Int=1)
    str = "ϱ"*subscriptify(i)*subscriptify(j)
    return AutSymbol(str, pow, RTransvect(i, j))
end

function lmul_autsymbol(i, j; pow::Int=1)
    str = "λ"*subscriptify(i)*subscriptify(j)
    return AutSymbol(str, pow, LTransvect(i, j))
end

function flip_autsymbol(i; pow::Int=1)
    pow = (2+pow%2)%2
    if pow == 0
       return id_autsymbol()
    else
        str = "ɛ"*subscriptify(i)
        return AutSymbol(str, pow, FlipAut(i))
    end
end

function perm_autsymbol(p::Generic.perm{Int8}; pow::Int=1)
    p = p^pow
    if p == parent(p)()
       return id_autsymbol()
    else
       str = "σ"*join([subscriptify(i) for i in p.d])
       return AutSymbol(str, 1, PermAut(p))
    end
end

function perm_autsymbol(a::Vector{T}) where T<:Integer
   G = PermutationGroup(Int8(length(a)))
   return perm_autsymbol(G(Vector{Int8}(a)))
end

domain(G::AutGroup)= NTuple{length(G.objectGroup.gens), FreeGroupElem}(gens(G.objectGroup))

###############################################################################
#
#   AutGroup / Automorphism constructors
#
###############################################################################

function AutGroup(G::FreeGroup; special=false)
   n = length(gens(G))
   n == 0 && return AutGroup{n}(G, AutSymbol[])
   S = AutSymbol[]

   indexing = [[i,j] for i in 1:n for j in 1:n if i≠j]

   rmuls = [rmul_autsymbol(i,j) for (i,j) in indexing]
   lmuls = [lmul_autsymbol(i,j) for (i,j) in indexing]

   append!(S, [rmuls; lmuls])

   if !special
      flips = [flip_autsymbol(i) for i in 1:n]
      syms = [perm_autsymbol(p) for p in elements(PermutationGroup(Int8(n)))][2:end]

      append!(S, [flips; syms])

   end
   return AutGroup{n}(G, S)
end

###############################################################################
#
#   Types call overloads
#
###############################################################################

function convert(::Type{Automorphism{N}}, s::AutSymbol) where N
    return Automorphism{N}(AutSymbol[s])
end

function (G::AutGroup{N})() where N
   id = Automorphism{N}(id_autsymbol())
   id.parent = G
   return id
end

function (G::AutGroup{N})(f::AutSymbol) where N
   g = Automorphism{N}([f])
   g.parent = G
   return g
end

function (G::AutGroup{N})(g::Automorphism{N}) where N
   g.parent = G
   return g
end

###############################################################################
#
#   Functional call overloads for evaluation of AutSymbol and Automorphism
#
###############################################################################

function (f::AutSymbol)(v::NTuple{N, T}) where {N, T}
   if f.pow == 0
      nothing
   else
      v = f.typ(v, f.pow)::NTuple{N, T}
   end
   return v
end

function (F::Automorphism{N})(v::NTuple{N, T}) where {N, T}
    for f in F.symbols
        v = f(v)::NTuple{N, T}
    end
    return v
end

###############################################################################
#
#   Comparison
#
###############################################################################

const HASHINGCONST = 0x314a5215305f3ec1 # more or less random

hash(s::AutSymbol, h::UInt) = hash(s.str, hash(s.pow, hash(:AutSymbol, h)))

function hash(g::Automorphism, h::UInt)
    if g.modified
        g.savedhash = hash(reduce!.(g(domain(parent(g)))),
            hash(typeof(g), hash(parent(g), HASHINGCONST)))
        g.modified = false
        end
    return xor(g.savedhash, h)
end

(==)(s::AutSymbol, t::AutSymbol) = s.str == t.str && s.pow == t.pow

function (==)(g::Automorphism{N}, h::Automorphism{N}) where N
    parent(g) == parent(h) || return false

    if !g.modified && !h.modified
        if g.savedhash != h.savedhash
            return false
        end
    end

    # expensive:
    g_im = reduce!.(g(domain(parent(g))))
    h_im = reduce!.(h(domain(parent(h))))
    # cheap:
    g.savedhash = hash(g_im,
        hash(typeof(g), hash(parent(g), HASHINGCONST)))
    g.modified = false
    h.savedhash = hash(h_im,
        hash(typeof(h), hash(parent(h), HASHINGCONST)))
    h.modified = false

    return g_im == h_im
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function change_pow(s::AutSymbol, n::Int)
    if n == 0
        return id_autsymbol()
    end
    symbol = s.typ
    if isa(symbol, FlipAut)
        return flip_autsymbol(symbol.i, pow=n)
    elseif isa(symbol, PermAut)
        return perm_autsymbol(symbol.perm, pow=n)
    elseif isa(symbol, RTransvect)
        return rmul_autsymbol(symbol.i, symbol.j, pow=n)
    elseif isa(symbol, LTransvect)
        return lmul_autsymbol(symbol.i, symbol.j, pow=n)
    elseif isa(symbol, Identity)
        return s
    else
        warn("Changing power of an unknown type of symbol! $s")
        return AutSymbol(s.str, n, s.typ)
    end
end

length(s::AutSymbol) = abs(s.pow)

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::AutGroup)
   print(io, "Automorphism Group of $(G.objectGroup)\n")
   print(io, "Generated by $(join(G.gens, ","))")
end

###############################################################################
#
#   Binary operators
#
###############################################################################

###############################################################################
#
#   Inversion
#
###############################################################################

inv(f::AutSymbol) = change_pow(f, -f.pow)

###############################################################################
#
#   Misc
#
###############################################################################

function getperm(s::AutSymbol)
    if s.pow != 1
        warn("Power for perm_symbol should be never 0!")
        return s.typ.perm^s.pow
    else
        return s.typ.perm
    end
end

function simplifyperms!(W::Automorphism{N}) where N
    reduced = true
    to_delete = Int[]
    for i in 1:length(W.symbols)-1
        if W.symbols[i].pow == 0
            continue
        elseif isa(W.symbols[i].typ, PermAut) && isa(W.symbols[i+1].typ, PermAut)
            reduced = false
            c = W.symbols[i]
            n = W.symbols[i+1]
            W.symbols[i+1] = perm_autsymbol(getperm(c)*getperm(n))
            push!(to_delete, i)
        end
    end
    deleteat!(W.symbols, to_delete)
    deleteids!(W)
    return reduced
end

function reduce!(W::Automorphism)
    if length(W) == 0
        return W
    elseif length(W.symbols) == 1
        deleteids!(W)
    else
        reduced = false
        while !reduced
            reduced = simplifyperms!(W) && freereduce!(W)
        end
    end

    W.modified = true

    return W
end
