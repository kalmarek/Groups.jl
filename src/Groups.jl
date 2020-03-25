module Groups

using AbstractAlgebra
import AbstractAlgebra: Group, GroupElem, Ring
import AbstractAlgebra: parent, parent_type, elem_type
import AbstractAlgebra: order, gens, matrix_repr

import Base: length, ==, hash, show, convert, eltype, iterate
import Base: inv, reduce, *, ^, power_by_squaring
import Base: findfirst, findnext
import Base: deepcopy_internal

export elements

using LinearAlgebra
using Markdown


Base.one(G::Generic.PermGroup) = Generic.Perm(G.n)
Base.one(r::NCRingElem) = one(parent(r))

###############################################################################
#
#   ParentType / ObjectType definition
#

abstract type AbstractFPGroup <: Group end

function Base.one(G::Gr) where Gr <: AbstractFPGroup
    El = elem_type(G)
    id = El(eltype(El)[])
    id.parent = G
    return id
end

elem_type(G::Gr) where Gr <:AbstractFPGroup = elem_type(Gr) # fallback definition

@doc doc"""
    ::GSymbol
> Abstract type which all group symbols of AbstractFPGroups should subtype. Each
> concrete subtype should implement fields:
> * `id` which is the `Symbol` representation/identification of a symbol
> * `pow` which is the (multiplicative) exponent of a symbol.

"""
abstract type GSymbol end

Base.iterate(s::GS, i=1) where GS<:GSymbol = i <= abs(s.pow) ? (GS(s.id, sign(s.pow)), i+1) : nothing
Base.length(s::GSymbol) = abs(s.pow)
Base.size(s::GSymbol) = (length(s), )
Base.eltype(s::GS) where GS<:GSymbol = GS
Base.isone(s::GSymbol) = iszero(s.pow)

change_pow(s::S, n::Integer) where S<:GSymbol = S(s.id, n)
Base.inv(s::GSymbol) = change_pow(s, -s.pow)

hash(s::S, h::UInt) where S<:GSymbol = hash(s.id, hash(s.pow, hash(S, h)))

abstract type GWord{T<:GSymbol} <: GroupElem end

# fallback definitions
Base.eltype(w::GW) where GW<:GWord = eltype(GW)
@doc doc"""
    W::GroupWord{T} <: GWord{T<:GSymbol} <:GroupElem
> Basic representation of element of a finitely presented group. `W.symbols`
> fieldname contains particular group symbols which multiplied constitute a
> group element, i.e. a word in generators.
> As reduction (inside group) of such word may be time consuming we provide
> `savedhash` and `modified` fields as well:
> hash (used e.g. in the `unique` function) is calculated by reducing the word,
> setting `modified` flag to `false` and computing the hash which is stored in
> `savedhash` field.
> whenever word `W` is changed `W.modified` is set to `false`;
> Future comparisons don't perform reduction (and use `savedhash`) as long as
> `modified` flag remains `false`.

"""
mutable struct GroupWord{T} <: GWord{T}
    symbols::Vector{T}
    modified::Bool
    savedhash::UInt
    parent::Group

    function GroupWord{T}(symbols::Vector{T}) where {T}
       return new{T}(symbols, true, zero(UInt))
    end
end

syllablelength(w::GWord) = length(w.symbols)
syllables(w::GWord) = w.symbols
ismodified(w::GWord) = w.modified
setmodified!(w::GWord) = (w.modified = true; w)
unsetmodified!(w::GWord) = (w.modified = false; w)
Base.one(w::GWord) = one(parent(w))


###############################################################################
#
#   Includes
#
###############################################################################

include("FreeGroup.jl")
include("FPGroups.jl")
include("AutGroup.jl")

include("DirectPower.jl")
include("WreathProducts.jl")

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent(w::GWord{T}) where {T<:GSymbol} = w.parent

###############################################################################
#
#   ParentType / ObjectType constructors
#
###############################################################################

GroupWord(s::T) where {T<:GSymbol} = GroupWord{T}(T[s])
GroupWord{T}(s::T) where {T<:GSymbol} = GroupWord{T}(T[s])
GroupWord(w::GroupWord{T}) where {T<:GSymbol} = w
convert(::Type{GroupWord{T}}, s::T) where {T<:GSymbol} = GroupWord{T}(T[s])

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function hash_internal(W::GWord)
    reduce!(W)
    return hash(syllables(W), hash(typeof(W), hash(parent(W))))
end

function hash(W::GWord, h::UInt)
    if ismodified(W)
        W.savedhash = hash_internal(W)
        unsetmodified!(W)
    end
    return xor(W.savedhash, h)
end

# WARNING: Due to specialised (constant) hash function of GWords this one is actually necessary!
function deepcopy_internal(W::T, dict::IdDict) where {T<:GWord}
    G = parent(W)
    return G(T(deepcopy(syllables(W))))
end

function freereduce!(::Type{Bool}, w::GWord)
    reduced = true
    for i in 1:syllablelength(w)-1
        s, ns = syllables(w)[i], syllables(w)[i+1]
        if isone(s)
            continue
        elseif s.id == ns.id
            reduced = false
            setmodified!(w)
            p1 = s.pow
            p2 = ns.pow

            syllables(w)[i+1] = change_pow(s, p1 + p2)
            syllables(w)[i] = change_pow(s, 0)
        end
    end
    filter!(!isone, syllables(w))
    return reduced
end

function freereduce!(w::GWord)
    reduced = false
    while !reduced
        reduced = freereduce!(Bool, w)
    end
    return w
end

reduce!(w::GWord) = freereduce!(w)

@doc doc"""
    reduce(w::GWord)
> performs reduction/simplification of a group element (word in generators).
> The default reduction is the free group reduction
> More specific procedures should be dispatched on `GWord`s type parameter.

"""
reduce(w::GWord) = reduce!(deepcopy(w))

@doc doc"""
    gens(G::AbstractFPGroups)
> returns vector of generators of `G`, as its elements.

"""
gens(G::AbstractFPGroup) = [G(g) for g in G.gens]

###############################################################################
#
#   String I/O
#
###############################################################################

@doc doc"""
    show(io::IO, W::GWord)
> The actual string produced by show depends on the eltype of `W.symbols`.

"""
function show(io::IO, W::GWord)
    if length(W) == 0
        print(io, "(id)")
    else
        join(io, [string(s) for s in W.symbols], "*")
    end
end

function show(io::IO, s::T) where {T<:GSymbol}
   if s.pow == 1
      print(io, string(s.id))
   else
      print(io, string((s.id))*"^$(s.pow)")
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

function (==)(W::T, Z::T) where T <: GWord
    parent(W) != parent(Z) && return false
    hash(W) != hash(Z) && return false
    return syllables(W) == syllables(Z)
end

function (==)(s::GSymbol, t::GSymbol)
    isone(s) && isone(t) && return true
    s.pow == t.pow && s.id == t.id && return true
    return false
end

###############################################################################
#
#   Binary operators
#
###############################################################################

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
        return one(parent(W))
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

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(W::T) where T<:GWord
    if length(W) == 0
        return W
    else
        G = parent(W)
        w = T([inv(s) for s in Iterators.reverse(syllables(W))])
        return G(w)
    end
end

###############################################################################
#
#   Replacement of symbols / sub-words
#
###############################################################################

issubsymbol(s::GSymbol, t::GSymbol) =
    s.id == t.id && (0 ≤ s.pow ≤ t.pow || 0 ≥ s.pow ≥ t.pow)

"""doc
Find the first linear index k>=i such that Z < W.symbols[k:k+length(Z)-1]
"""
function findnext(W::GWord, Z::GWord, i::Int)
    n = length(Z.symbols)
   if n == 0
      return 0
   elseif n == 1
      for idx in i:lastindex(W.symbols)
         if issubsymbol(Z.symbols[1], W.symbols[idx])
            return idx
         end
      end
      return 0
   else
      for idx in i:lastindex(W.symbols) - n + 1
         foundfirst = issubsymbol(Z.symbols[1], W.symbols[idx])
         lastmatch = issubsymbol(Z.symbols[end], W.symbols[idx+n-1])
         if foundfirst && lastmatch
            # middles match:
            if view(Z.symbols, 2:n-1) == view(W.symbols, idx+1:idx+n-2)
               return idx
            end
         end
      end
   end
   return 0
end

findfirst(W::GWord, Z::GWord) = findnext(W, Z, 1)

function replace!(W::GWord, index, toreplace::GWord, replacement::GWord; check=true)
   n = length(toreplace.symbols)
   if n == 0
      return reduce!(W)

   elseif n == 1
      if check
         @assert issubsymbol(toreplace.symbols[1], W.symbols[index])
      end

      first = change_pow(W.symbols[index],
         W.symbols[index].pow - toreplace.symbols[1].pow)
      last = change_pow(W.symbols[index], 0)

   else
      if check
         @assert issubsymbol(toreplace.symbols[1], W.symbols[index])
         @assert W.symbols[index+1:index+n-2] == toreplace.symbols[2:end-1]
         @assert issubsymbol(toreplace.symbols[end], W.symbols[index+n-1])
      end

      first = change_pow(W.symbols[index],
      W.symbols[index].pow - toreplace.symbols[1].pow)
      last = change_pow(W.symbols[index+n-1],
      W.symbols[index+n-1].pow - toreplace.symbols[end].pow)
   end

   replacement = first * replacement * last
   splice!(W.symbols, index:index+n-1, replacement.symbols)
   return reduce!(W)
end

function replace(W::GWord, index, toreplace::GWord, replacement::GWord)
    replace!(deepcopy(W), index, toreplace, replacement)
end

function replace_all!(W::T,subst_dict::Dict{T,T}) where {T<:GWord}
    modified = false
    for toreplace in reverse!(sort!(collect(keys(subst_dict)), by=length))
        replacement = subst_dict[toreplace]
        i = findfirst(W, toreplace)
        while i ≠ 0
            modified = true
            replace!(W,i,toreplace, replacement)
            i = findnext(W, toreplace, i)
        end
    end
    return modified
end

function replace_all(W::T, subst_dict::Dict{T,T}) where {T<:GWord}
   W = deepcopy(W)
   replace_all!(W, subst_dict)
   return W
end

###############################################################################
#
#   Misc
#
###############################################################################

function generate_balls(S::AbstractVector{T}, Id::T=one(parent(first(S)));
        radius=2, op=*) where T<:GroupElem
    sizes = Int[]
    B = [Id]
    for i in 1:radius
        BB = [op(i,j) for (i,j) in Base.product(B,S)]
        B = unique([B; vec(BB)])
        push!(sizes, length(B))
    end
    return B, sizes
end

function generate_balls(S::AbstractVector{T}, Id::T=one(parent(first(S)));
        radius=2, op=*) where {T<:NCRingElem}
    sizes = Int[]
    B = [Id]
    for i in 1:radius
        BB = [op(i,j) for (i,j) in Base.product(B,S)]
        B = unique([B; vec(BB)])
        push!(sizes, length(B))
    end
    return B, sizes
end

########### iteration for GFField


length(F::AbstractAlgebra.GFField) = order(F)

function iterate(F::AbstractAlgebra.GFField, s=0)
   if s >= order(F)
      return nothing
   else
      return F(s), s+1
   end
end

eltype(::Type{AbstractAlgebra.GFField{I}}) where I = AbstractAlgebra.gfelem{I}

end # of module Groups
