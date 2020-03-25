###############################################################################
#
#   AutSymbol/ AutGroup / Automorphism
#
###############################################################################

struct RTransvect
    i::Int8
    j::Int8
end

struct LTransvect
    i::Int8
    j::Int8
end

struct FlipAut
    i::Int8
end

struct PermAut
    perm::Generic.Perm{Int8}
end

struct Identity end

struct AutSymbol <: GSymbol
    id::Symbol
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

    function Automorphism{N}(f::Vector{AutSymbol}) where {N}
        return new{N}(f, true, zero(UInt))
    end
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

function (ϱ::RTransvect)(v, pow::Integer=1)
    append!(v[ϱ.i], v[ϱ.j]^pow)
    freereduce!(v[ϱ.i])
    return v
end

function (λ::LTransvect)(v, pow::Integer=1)
    prepend!(v[λ.i], v[λ.j]^pow)
    freereduce!(v[λ.i])
    return v
end

function (σ::PermAut)(v, pow::Integer=1)
    w = deepcopy(v)
    s = (σ.perm^pow).d
    @inbounds for k in eachindex(v)
        v[k].symbols = w[s[k]].symbols
    end
    return v
end

function (ɛ::FlipAut)(v, pow::Integer=1)
    @inbounds if isodd(pow)
        v[ɛ.i].symbols = inv(v[ɛ.i]).symbols
    end
    return v
end

(::Identity)(v, pow::Integer=1) = v

# taken from ValidatedNumerics, under under the MIT "Expat" License:
# https://github.com/JuliaIntervals/ValidatedNumerics.jl/blob/master/LICENSE.md
function subscriptify(n::Integer)
    subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
    @assert 0 <= n <= 9
    return Char(subscript_0 + n)
    # return [Char(subscript_0 + i) for i in reverse(digits(n))])
end

function id_autsymbol()
   return AutSymbol(Symbol("(id)"), 0, Identity())
end

function transvection_R(i::Integer, j::Integer, pow::Integer=1)
    id = Symbol("ϱ", subscriptify(i), subscriptify(j))
    return AutSymbol(id, pow, RTransvect(i, j))
end

function transvection_L(i::Integer, j::Integer, pow::Integer=1)
    id = Symbol("λ", subscriptify(i), subscriptify(j))
    return AutSymbol(id, pow, LTransvect(i, j))
end

function flip(i::Integer, pow::Integer=1)
    iseven(pow) && return id_autsymbol()
    id = Symbol("ɛ", subscriptify(i))
    return AutSymbol(id, 1, FlipAut(i))
end

function AutSymbol(p::Generic.Perm, pow::Integer=1)
    if pow != 1
        p = p^pow
    end

    if any(p.d[i] != i for i in eachindex(p.d))
        id = Symbol("σ", "₍", [subscriptify(i) for i in p.d]..., "₎")
        return AutSymbol(id, 1, PermAut(p))
    end
    return id_autsymbol()
end

function AutSymbol(a::Vector{<:Integer}, pow=1)
   return AutSymbol(Generic.Perm(convert(Vector{Int8}, a)), pow)
end

ϱ(i::Integer, j::Integer, pow::Integer=1) = transvection_R(i, j, pow=pow)
λ(i::Integer, j::Integer, pow::Integer=1) = transvection_L(i, j, pow=pow)
ε(i::Integer, pow::Integer=1) = flip(i, pow=pow)
σ(v, pow=1) = AutSymbol(v, pow=pow)

function domain(G::AutGroup{N}) where N
    F = G.objectGroup
    gg = gens(F)
    return ntuple(i->gg[i], N)
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

   indexing = [[i,j] for i in 1:n for j in 1:n if i≠j]

   rmuls = [transvection_R(i,j) for (i,j) in indexing]
   lmuls = [transvection_L(i,j) for (i,j) in indexing]

   append!(S, [rmuls; lmuls])

   if !special
      flips = [flip(i) for i in 1:n]
      syms = [AutSymbol(p) for p in PermutationGroup(Int8(n))][2:end]

      append!(S, [flips; syms])
   end
   return AutGroup{n}(G, S)
end

Aut(G::Group) = AutGroup(G)
SAut(G::Group) = AutGroup(G, special=true)

###############################################################################
#
#   Types call overloads
#
###############################################################################

Automorphism{N}(s::AutSymbol) where N = Automorphism{N}(AutSymbol[s])

function Base.one(G::AutGroup{N}) where N
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

(s::AutSymbol)(v::NTuple{N, T}) where {N, T} = s.fn(v, s.pow)::NTuple{N, T}

function (f::Automorphism{N})(v::NTuple{N, T}) where {N, T}
    for s in syllables(f)
        v = s(v)::NTuple{N, T}
        # if iszero(i % 3)
            # freereduce!.(v)
        # end
    end
    # return freereduce!.(v)
    return v
end

evaluate(f::Automorphism) = f(domain(parent(f)))

###############################################################################
#
#   Comparison
#
###############################################################################

const HASHINGCONST = 0x7d28276b01874b19 # hash(Automorphism)

hash(s::AutSymbol, h::UInt) = hash(s.id, hash(s.pow, hash(AutSymbol, h)))

function hash_internal(g::Automorphism, images = freereduce!.(evaluate(g)),
    h::UInt = HASHINGCONST)
    return hash(images, hash(parent(g), h))
end

function compute_images(g::Automorphism)
    images = reduce!.(evaluate(g))
    g.savedhash = hash_internal(g, images)
    unsetmodified!(g)
    return images
end

function (==)(g::Automorphism{N}, h::Automorphism{N}) where N
    img_c, imh_c = false, false

    if ismodified(g)
        img = compute_images(g)
        img_c = true
    end

    if ismodified(h)
        imh = compute_images(h)
        imh_c = true
    end

    @assert !ismodified(g) && !ismodified(h)
    # cheap
    hash(g) != hash(h) && return false # hashes differ, so images must have differed as well
    # equal elements, or possibly hash conflict
    if !img_c
        img = compute_images(g)
    end
    if !imh_c
        imh = compute_images(h)
    end
    return img == imh
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
        return flip(symbol.i, n)
    elseif symbol isa PermAut
        return AutSymbol(symbol.perm, n)
    elseif symbol isa RTransvect
        return transvection_R(symbol.i, symbol.j, n)
    elseif symbol isa LTransvect
        return transvection_L(symbol.i, symbol.j, n)
    elseif symbol isa Identity
        return s
    else
        throw(DomainError("Unknown type of AutSymbol: $s"))
    end
end

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

###############################################################################
#
#   Misc
#
###############################################################################

getperm(s::AutSymbol) = s.fn.perm^s.pow

function simplifyperms!(::Type{Bool}, w::Automorphism{N}) where N
    reduced = true
    for i in 1:syllablelength(w)-1
        s, ns = syllables(w)[i], syllables(w)[i+1]
        if isone(s)
            continue
        elseif s.fn isa PermAut && ns.fn isa PermAut
            reduced = false
            setmodified!(w)
            syllables(w)[i+1] = AutSymbol(getperm(s)*getperm(ns))
            syllables(w)[i] = change_pow(s, 0)
        end
    end
    filter!(!isone, syllables(w))
    return reduced
end

function reduce!(w::Automorphism)
    reduced = false
    while !reduced
        reduced = simplifyperms!(Bool, w) && freereduce!(Bool, w)
    end
    return W
end

function linear_repr(A::Automorphism{N}, hom=matrix_repr) where N
    return reduce(*, linear_repr.(A.symbols, N, hom), init=hom(Identity(),N,1))
end

linear_repr(a::AutSymbol, n::Int, hom) = hom(a.fn, n, a.pow)

function matrix_repr(a::Union{RTransvect, LTransvect}, n::Int, pow)
    x = Matrix{Int}(I, n, n)
    x[a.i,a.j] = pow
    return x
end

function matrix_repr(a::FlipAut, n::Int, pow)
    x = Matrix{Int}(I, n, n)
    x[a.i,a.i] = -1^pow
    return x
end

matrix_repr(a::PermAut, n::Int, pow) = Matrix{Int}(I, n, n)[(a.perm^pow).d, :]

matrix_repr(a::Identity, n::Int, pow) = Matrix{Int}(I, n, n)
