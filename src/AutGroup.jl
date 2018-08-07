###############################################################################
#
#   AutSymbol/ AutGroup / Automorphism
#
###############################################################################

struct RTransvect{I<:Integer}
    i::I
    j::I
end

struct LTransvect{I<:Integer}
    i::I
    j::I
end

struct FlipAut{I<:Integer}
    i::I
end

struct PermAut{I<:Integer}
    perm::Generic.perm{I}
end

struct Identity end

struct AutSymbol <: GSymbol
   str::String
   pow::Int8
   fn::Union{LTransvect, RTransvect, PermAut, FlipAut, Identity}
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

    Automorphism{N}(f::Vector{AutSymbol}) where N = new(f, true, zero(UInt))

end

export Automorphism, AutGroup, Aut, SAut

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

function (ϱ::RTransvect{I})(v, pow::Integer=one(I)) where I
    @inbounds Groups.r_multiply!(v[ϱ.i], (v[ϱ.j]^pow).symbols, reduced=false)
    return v
end

function (λ::LTransvect{I})(v, pow::Integer=one(I)) where I
    @inbounds Groups.l_multiply!(v[λ.i], (v[λ.j]^pow).symbols, reduced=false)
    return v
end

function (σ::PermAut{I})(v, pow::Integer=one(I)) where I
   w = deepcopy(v)
   s = (σ.perm^pow).d
   @inbounds for k in eachindex(v)
       v[k].symbols = w[s[k]].symbols
   end
   return v
end

function (ɛ::FlipAut{I})(v, pow::Integer=one(I)) where I
   @inbounds if isodd(pow)
       v[ɛ.i].symbols = inv(v[ɛ.i]).symbols
   end
   return v
end

(::Identity)(v, pow::Integer=zero(Int8)) = v

# taken from ValidatedNumerics, under under the MIT "Expat" License:
# https://github.com/JuliaIntervals/ValidatedNumerics.jl/blob/master/LICENSE.md
function subscriptify(n::Integer)
    subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
    return join([Char(subscript_0 + i) for i in reverse(digits(n))])
end

function id_autsymbol()
   return AutSymbol("(id)", 0, Identity())
end

function rmul_autsymbol(i::I, j::I; pow::Integer=one(I)) where I<:Integer
    str = "ϱ"*subscriptify(i)*subscriptify(j)
    return AutSymbol(str, I(pow), RTransvect(i, j))
end

function lmul_autsymbol(i::I, j::I; pow::Integer=one(I)) where I<:Integer
    str = "λ"*subscriptify(i)*subscriptify(j)
    return AutSymbol(str, I(pow), LTransvect(i, j))
end

function flip_autsymbol(i::I; pow::Integer=one(I)) where I<:Integer
    pow = I((2+pow%2)%2)
    if pow == zero(I)
       return id_autsymbol()
    else
        str = "ɛ"*subscriptify(i)
        return AutSymbol(str, I(pow), FlipAut(i))
    end
end

function perm_autsymbol(p::Generic.perm{I}; pow::Integer=one(I)) where I<:Integer
    p = p^pow
    for i in eachindex(p.d)
        if p.d[i] != i
            str = "σ"*join([subscriptify(i) for i in p.d])
            return AutSymbol(str, one(I), PermAut(p))
        end
    end
    return id_autsymbol()
end

function perm_autsymbol(a::Vector{T}) where T<:Integer
   return perm_autsymbol(perm(Vector{Int8}(a), false))
end

function domain(G::AutGroup{N}) where N
    F = G.objectGroup
    gg = gens(F)
    return ntuple(i->gg[i], Val{N})
end

###############################################################################
#
#   AutGroup / Automorphism constructors
#
###############################################################################

function AutGroup(G::FreeGroup; special=false)
   S = AutSymbol[]
   n = length(gens(G))
   n == 0 && return AutGroup{n}(G, S)

   n = convert(Int8, n)

   indexing = [[i,j] for i in Int8(1):n for j in Int8(1):n if i≠j]

   rmuls = [rmul_autsymbol(i,j) for (i,j) in indexing]
   lmuls = [lmul_autsymbol(i,j) for (i,j) in indexing]

   append!(S, [rmuls; lmuls])

   if !special
      flips = [flip_autsymbol(i) for i in 1:n]
      syms = [perm_autsymbol(p) for p in elements(PermutationGroup(n))][2:end]

      append!(S, [flips; syms])

   end
   return AutGroup{Int64(n)}(G, S)
end

Aut(G::Group) = AutGroup(G)
SAut(G::Group) = AutGroup(G, special=true)

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
   if f.pow != 0
      v = f.fn(v, f.pow)::NTuple{N, T}
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

const HASHINGCONST = 0x7d28276b01874b19 # hash(Automorphism)

hash(s::AutSymbol, h::UInt) = hash(s.str, hash(s.pow, hash(:AutSymbol, h)))

function hash(g::Automorphism, h::UInt)
    if g.modified
        g_im = reduce!.(g(domain(parent(g))))
        g.savedhash = hash(g_im, hash(typeof(g), hash(parent(g), HASHINGCONST)))
        g.modified = false
        end
    return xor(g.savedhash, h)
end

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
    g.savedhash = hash(g_im, hash(typeof(g), hash(parent(g), HASHINGCONST)))
    g.modified = false
    h.savedhash = hash(h_im, hash(typeof(h), hash(parent(h), HASHINGCONST)))
    h.modified = false

    return g_im == h_im
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function change_pow(s::AutSymbol, n::Integer)
    if n == zero(n)
        return id_autsymbol()
    end
    symbol = s.fn
    if symbol isa FlipAut
        return flip_autsymbol(symbol.i, pow=n)
    elseif symbol isa PermAut
        return perm_autsymbol(symbol.perm, pow=n)
    elseif symbol isa RTransvect
        return rmul_autsymbol(symbol.i, symbol.j, pow=n)
    elseif symbol isa LTransvect
        return lmul_autsymbol(symbol.i, symbol.j, pow=n)
    elseif symbol isa Identity
        return s
    else
        warn("Changing power of an unknown type of symbol! $s")
        return AutSymbol(s.str, n, s.fn)
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
        return s.fn.perm^s.pow
    else
        return s.fn.perm
    end
end

function simplifyperms!(W::Automorphism{N}) where N
    reduced = true
    to_delete = Int[]
    for i in 1:length(W.symbols)-1
        if W.symbols[i].pow == 0
            continue
        elseif W.symbols[i].fn isa PermAut && W.symbols[i+1].fn isa PermAut
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

function linear_repr(A::Automorphism{N}, hom=matrix_repr) where N
    return reduce(*, hom(Identity(), N, 1), linear_repr.(A.symbols, N, hom))
end

linear_repr(a::AutSymbol, n::Int, hom) = hom(a.fn, n, a.pow)

function matrix_repr(a::Union{RTransvect, LTransvect}, n::Int, pow)
    x = eye(n)
    x[a.i,a.j] = pow
    return x
end

function matrix_repr(a::FlipAut, n::Int, pow)
    x = eye(n)
    x[a.i,a.i] = -1^pow
    return x
end

matrix_repr(a::PermAut, n::Int, pow) = eye(n)[(a.perm^pow).d, :]

matrix_repr(a::Identity, n::Int, pow) = eye(n)
