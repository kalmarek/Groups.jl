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

Base.one(G::Generic.PermGroup) = G(collect(1:G.n), false)

###############################################################################
#
#   ParentType / ObjectType definition
#
###############################################################################

@doc doc"""
    ::GSymbol
> Abstract type which all group symbols of AbstractFPGroups should subtype. Each
> concrete subtype should implement fields:
> * `id` which is the `Symbol` representation/identification of a symbol
> * `pow` which is the (multiplicative) exponent of a symbol.

"""
abstract type GSymbol end

abstract type GWord{T<:GSymbol} <:GroupElem end

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
   savedhash::UInt
   modified::Bool
   parent::Group

   function GroupWord{T}(symbols::Vector{T}) where {T}
      return new{T}(symbols, hash(symbols), true)
   end
end

abstract type AbstractFPGroup <: Group end

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

function hash(W::GWord, h::UInt)
    W.modified && reduce!(W)
    return xor(W.savedhash, h)
end

# WARNING: Due to specialised (constant) hash function of GWords this one is actually necessary!
function deepcopy_internal(W::T, dict::IdDict) where {T<:GWord}
    G = parent(W)
    return G(T(deepcopy(W.symbols)))
end

length(W::GWord) = sum([length(s) for s in W.symbols])

function deleteids!(W::GWord)
    to_delete = Int[]
    for i in 1:length(W.symbols)
        if W.symbols[i].pow == 0
           push!(to_delete, i)
        end
    end
    deleteat!(W.symbols, to_delete)
end

function freereduce!(W::GWord)
    reduced = true
    for i in 1:length(W.symbols) - 1
        if W.symbols[i].pow == 0
            continue
        elseif W.symbols[i].id == W.symbols[i+1].id
            reduced = false
            p1 = W.symbols[i].pow
            p2 = W.symbols[i+1].pow

            W.symbols[i+1] = change_pow(W.symbols[i], p1 + p2)
            W.symbols[i] = change_pow(W.symbols[i], 0)
        end
    end
    deleteids!(W)
    return reduced
end

function reduce!(W::GWord)
    if length(W) < 2
        deleteids!(W)
    else
        reduced = false
        while !reduced
            reduced = freereduce!(W)
        end
    end

    W.savedhash = hash(W.symbols, hash(typeof(W), hash(parent(W), zero(UInt))))
    W.modified = false

    return W
end

@doc doc"""
    reduce(W::GWord)
> performs reduction/simplification of a group element (word in generators).
> The default reduction is the free group reduction, i.e. consists of
> multiplying adjacent symbols with the same `id` identifier and deleting the
> identity elements from `W.symbols`.
> More specific procedures should be dispatched on `GWord`s type parameter.

"""
reduce(W::GWord) = reduce!(deepcopy(W))

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

function (==)(W::GWord, Z::GWord)
    parent(W) == parent(Z) || return false

    W.modified && reduce!(W)
    Z.modified && reduce!(Z)

    if W.savedhash != Z.savedhash
        return false
    end

    return W.symbols == Z.symbols
end

function (==)(s::GSymbol, t::GSymbol)
   s.pow == t.pow || return false
   s.pow ==  0 && return true
   s.id == t.id || return false
   return true
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function AbstractAlgebra.mul!(out::GWord, x::GWord, y::GWord; reduced::Bool=true)
    resize!(out.symbols, length(x.symbols)+length(y.symbols))
    for i in eachindex(x.symbols)
        out.symbols[i] = x.symbols[i]
    end
    for i in eachindex(y.symbols)
        out.symbols[length(x.symbols)+i] = y.symbols[i]
    end
    if reduced
        reduce!(out)
    end
    return out
end

function r_multiply!(W::GWord, x; reduced::Bool=true)
    if length(x) > 0
        append!(W.symbols, x)
    end
    if reduced
        reduce!(W)
    end
    return W
end

function l_multiply!(W::GWord, x; reduced::Bool=true)
    if length(x) > 0
        prepend!(W.symbols, x)
    end
    if reduced
        reduce!(W)
    end
    return W
end

r_multiply(W::GWord, x; reduced=true) =
    r_multiply!(deepcopy(W),x, reduced=reduced)
l_multiply(W::GWord, x; reduced=true) =
    l_multiply!(deepcopy(W),x, reduced=reduced)

(*)(W::GWord, Z::GWord) = r_multiply(W, Z.symbols)
(*)(W::GWord, s::GSymbol) = r_multiply(W, [s])
(*)(s::GSymbol, W::GWord) = l_multiply(W, [s])

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
        r_multiply!(W, W.symbols)
    end
    Z = deepcopy(W)
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) >= 0
            r_multiply!(W, W.symbols)
        end
        r_multiply!(Z, W.symbols)
    end
    return Z
end

(^)(x::GWord, n::Integer) = power_by_squaring(x,n)

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(W::T) where {T<:GWord}
    if length(W) == 0
        return W
    else
        G = parent(W)
        w = T(reverse([inv(s) for s in W.symbols]))
        w.modified = true
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
