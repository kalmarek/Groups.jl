export Automorphism, AutGroup, Aut, SAut

###############################################################################
#
#   AutSymbol/ AutGroup / Automorphism
#

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

# taken from ValidatedNumerics, under under the MIT "Expat" License:
# https://github.com/JuliaIntervals/ValidatedNumerics.jl/blob/master/LICENSE.md
function subscriptify(n::Integer)
    subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
    return join([Char(subscript_0 + i) for i in reverse(digits(n))], "")
end

function id_autsymbol()
   return AutSymbol(Symbol("(id)"), 0, Identity())
end

function transvection_R(i::Integer, j::Integer, pow::Integer=1)
    if 0 < i < 10 && 0 < j < 10
        id = Symbol(:ϱ, subscriptify(i), subscriptify(j))
    else
        id = Symbol(:ϱ, subscriptify(i), "." ,subscriptify(j))
    end
    return AutSymbol(id, pow, RTransvect(i, j))
end

function transvection_L(i::Integer, j::Integer, pow::Integer=1)
    if 0 < i < 10 && 0 < j < 10
        id = Symbol(:λ, subscriptify(i), subscriptify(j))
    else
        id = Symbol(:λ, subscriptify(i), "." ,subscriptify(j))
    end
    return AutSymbol(id, pow, LTransvect(i, j))
end

function flip(i::Integer, pow::Integer=1)
    iseven(pow) && return id_autsymbol()
    id = Symbol(:ɛ, subscriptify(i))
    return AutSymbol(id, 1, FlipAut(i))
end

function AutSymbol(p::Generic.Perm, pow::Integer=1)
    if pow != 1
        p = p^pow
    end

    if any(p.d[i] != i for i in eachindex(p.d))
        id = Symbol(:σ, "₍", join([subscriptify(i) for i in p.d],""), "₎")
        return AutSymbol(id, 1, PermAut(p))
    end
    return id_autsymbol()
end

ϱ(i::Integer, j::Integer, pow::Integer=1) = transvection_R(i, j, pow)
λ(i::Integer, j::Integer, pow::Integer=1) = transvection_L(i, j, pow)
ε(i::Integer, pow::Integer=1) = flip(i, pow)
σ(v::Generic.Perm, pow::Integer=1) = AutSymbol(v, pow)

function change_pow(s::AutSymbol, n::Integer)
    iszero(n) && id_autsymbol()

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
        return id_autsymbol()
    else
        throw(DomainError("Unknown type of AutSymbol: $s"))
    end
end

###############################################################################
#
#   AutGroup / Automorphism
#

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

elem_type(::Type{AutGroup{N}}) where N = Automorphism{N}
parent_type(::Type{Automorphism{N}}) where N = AutGroup{N}

function AutGroup(G::FreeGroup; special=false)
   S = AutSymbol[]
   n = length(gens(G))
   n == 0 && return AutGroup{n}(G, S)

   indexing = [[i,j] for i in 1:n for j in 1:n if i≠j]

   rmuls = [ϱ(i,j) for (i,j) in indexing]
   lmuls = [λ(i,j) for (i,j) in indexing]

   append!(S, [rmuls; lmuls])

   if !special
      flips = [ε(i) for i in 1:n]
      syms = [σ(p) for p in SymmetricGroup(Int8(n))][2:end]

      append!(S, [flips; syms])
   end
   return AutGroup{n}(G, S)
end

Aut(G::Group) = AutGroup(G)
SAut(G::Group) = AutGroup(G, special=true)

Automorphism{N}(s::AutSymbol) where N = Automorphism{N}(AutSymbol[s])

function (G::AutGroup{N})(f::AutSymbol) where N
   g = Automorphism{N}([f])
   setparent!(g, G)
   return g
end

(G::AutGroup{N})(g::Automorphism{N}) where N = (setparent!(g, G); g)

###############################################################################
#
#   AutSymbol defining functions && evaluation
#   NOTE: all automorphisms operate on a tuple of FreeWords INPLACE!
#

function (ϱ::RTransvect)(v, pow::Integer=1)
    rmul!(v[ϱ.i], v[ϱ.j]^pow)
    return v
end

function (λ::LTransvect)(v, pow::Integer=1)
    lmul!(v[λ.i], v[λ.j]^pow)
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

###############################################################################
#
#   Functional call overloads for evaluation of AutSymbol and Automorphism
#

(s::AutSymbol)(v::NTuple{N, T}) where {N, T} = s.fn(v, s.pow)::NTuple{N, T}

function (f::Automorphism{N})(v::NTuple{N, T}) where {N, T}
    for s in syllables(f)
        v = s(v)::NTuple{N, T}
    end
    return v
end

function domain(G::AutGroup{N}) where N
    F = G.objectGroup
    return ntuple(i->F(F.gens[i]), N)
end

evaluate(f::Automorphism) = f(domain(parent(f)))

###############################################################################
#
#   hashing && equality
#

function hash_internal(
    g::Automorphism,
    h::UInt = 0x7d28276b01874b19; # hash(Automorphism)
    # alternatively: 0xcbf29ce484222325 from FNV-1a algorithm
    images = compute_images(g),
    prime = 0x00000100000001b3, # prime from FNV-1a algorithm
)
    return foldl((h,x) -> hash(x, h)*prime, images, init = hash(parent(g), h))
end

function compute_images(g::Automorphism)
    images = evaluate(g)
    for im in images
        reduce!(im)
    end
    return images
end

function (==)(g::Automorphism{N}, h::Automorphism{N}) where N
    syllables(g) == syllables(h) && return true
    img_computed, imh_computed = false, false

    if ismodified(g)
        img = compute_images(g) # sets modified bit
        hash(g, images=img)
        img_computed = true
    end

    if ismodified(h)
        imh = compute_images(h) # sets modified bit
        hash(h, images=imh)
        imh_computed = true
    end

    @assert !ismodified(g) && !ismodified(h)
    # cheap
    # if hashes differ, images must have differed as well
    hash(g) != hash(h) && return false

    # hashes equal, hence either equal elements, or a hash conflict
    begin
        if !img_computed
            img_task = Threads.@spawn img = compute_images(g)
            # img = compute_images(g)
        end
        if !imh_computed
            imh_task = Threads.@spawn imh = compute_images(h)
            # imh = compute_images(h)
        end
        !img_computed && fetch(img_task)
        !imh_computed && fetch(imh_task)
    end

    img != imh && @warn "hash collision in == :" g h
    return img == imh
end

###############################################################################
#
#   String I/O
#

function show(io::IO, G::AutGroup)
   print(io, "Automorphism Group of $(G.objectGroup)\n")
   print(io, "Generated by $(gens(G))")
end

###############################################################################
#
#   Reduction
#

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
    return w
end

###############################################################################
#
#   Abelianization (natural Representation to GL(N,Z))
#

abelianize(A::Automorphism{N}) where N = image(A, abelianize; n=N)

# homomorphism definition
abelianize(; n::Integer=1) = Matrix{Int}(I, n, n)
abelianize(a::AutSymbol; n::Int=1) = abelianize(a.fn, n, a.pow)

function abelianize(a::Union{RTransvect, LTransvect}, n::Int, pow)
    x = Matrix{Int}(I, n, n)
    x[a.i,a.j] = pow
    return x
end

function abelianize(a::FlipAut, n::Int, pow)
    x = Matrix{Int}(I, n, n)
    x[a.i,a.i] = -1
    return x
end

abelianize(a::PermAut, n::Integer, pow) = Matrix{Int}(I, n, n)[(a.perm^pow).d, :]
abelianize(a::Identity, n::Integer, pow) = abelianize(;n=n)
