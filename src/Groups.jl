__precompile__()
module Groups

using Nemo
import Nemo: Group, GroupElem, Ring
import Nemo: parent, parent_type, elem_type
import Nemo: elements, order, gens

import Base: length, ==, hash, show, convert
import Base: inv, reduce, *, ^
import Base: findfirst, findnext
import Base: deepcopy_internal

###############################################################################
#
#   ParentType / ObjectType definition
#
###############################################################################

doc"""
    ::GSymbol
> Abstract type which all group symbols of AbstractFPGroups should subtype. Each
> concrete subtype should implement fields:
> * `str` which is the string representation/identification of a symbol
> * `pow` which is the (multiplicative) exponent of a symbol.

"""
abstract type GSymbol end

doc"""
    W::GWord{T<:GSymbol} <:GroupElem
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

abstract type GWord{T<:GSymbol} <:GroupElem end

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
#   Type and parent object methods
#
###############################################################################

parent{T<:GSymbol}(w::GWord{T}) = w.parent

###############################################################################
#
#   ParentType / ObjectType constructors
#
###############################################################################

GroupWord(s::T) where {T<:GSymbol} = GroupWord{T}(T[s])
convert(::Type{GroupWord{T}}, s::T) where {T<:GSymbol} = GroupWord{T}(T[s])

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function hash(W::GWord, h::UInt)
    W.modified && reduce!(W)
    res = xor(W.savedhash, h)
    return res
end

function deepcopy_internal(W::T, dict::ObjectIdDict) where {T<:GWord}
    G = parent(W)
    return G(T(deepcopy(W.symbols)))
end

isone(s::GSymbol) = s.pow == 0

length(W::GWord) = sum([length(s) for s in W.symbols])

function free_reduce!(W::GWord)
    reduced = true
    for i in 1:length(W.symbols) - 1
        if W.symbols[i].str == W.symbols[i+1].str
            reduced = false
            p1 = W.symbols[i].pow
            p2 = W.symbols[i+1].pow
            W.symbols[i+1] = change_pow(W.symbols[i], p1 + p2)
            W.symbols[i] = change_pow(W.symbols[i], 0)
        end
    end
    deleteat!(W.symbols, find(x -> x.pow == 0, W.symbols))
    return reduced
end

function reduce!(W::GWord)
    if length(W) < 2
        deleteat!(W.symbols, find(x -> x.pow == 0, W.symbols))
    else
        reduced = false
        while !reduced
            reduced = free_reduce!(W)
        end
    end

    W.modified = false
    W.savedhash = hash(W.symbols, hash(typeof(W)))
    return W
end

doc"""
    reduce(W::GWord)
> performs reduction/simplification of a group element (word in generators).
> The default reduction is the free group reduction, i.e. consists of
> multiplying adjacent symbols with the same `str` identifier and deleting the
> identity elements from `W.symbols`.
> More specific procedures should be dispatched on `GWord`s type parameter.

"""
reduce(W::GWord) = reduce!(deepcopy(W))

doc"""
    gens(G::AbstractFPGroups)
> returns vector of generators of `G`, as its elements.

"""
gens(G::AbstractFPGroup) = [G(g) for g in G.gens]

###############################################################################
#
#   String I/O
#
###############################################################################

doc"""
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
   if isone(s)
      print(io, "(id)")
   elseif s.pow == 1
      print(io, s.str)
   else
      print(io, (s.str)*"^$(s.pow)")
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

function (==)(W::GWord, Z::GWord)
    parent(W) == parent(Z) || return false
    W.modified && reduce!(W) # reduce clears the flag and calculates savedhash
    Z.modified && reduce!(Z)
    return W.savedhash == Z.savedhash && W.symbols == Z.symbols
end

###############################################################################
#
#   Binary operators
#
###############################################################################

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
        return parent(W)()
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

is_subsymbol(s::GSymbol, t::GSymbol) =
    s.str == t.str && (0 ≤ s.pow ≤ t.pow || 0 ≥ s.pow ≥ t.pow)

function findfirst(W::GWord, Z::GWord)
    n = length(Z.symbols)

   if n == 0
      return 0
   elseif n == 1
      for (i,s) in enumerate(W.symbols)
         if is_subsymbol(Z.symbols[1], s)
            return i
         end
      end
      return 0
   else
      for (idx,a) in enumerate(W.symbols)
         if idx + n - 1 > length(W.symbols)
            break
         end
         foundfirst = is_subsymbol(Z.symbols[1], a)
         if foundfirst
            middlematch = W.symbols[idx+1:idx+n-2] == Z.symbols[2:end-1]
            lastmatch = is_subsymbol(Z.symbols[end], W.symbols[idx+n-1])
            if middlematch && lastmatch
               return idx
            end
         end
      end
   end

   return 0
end

function findnext(W::GWord, Z::GWord, i::Integer)
    t = findfirst(GWord{eltype(W.symbols)}(W.symbols[i:end]), Z)
    if t > 0
        return t+i-1
    else
        return 0
    end
end

function replace!(W::GWord, index, toreplace::GWord, replacement::GWord; check=true)
   n = length(toreplace.symbols)
   if n == 0
      return reduce!(W)

   elseif n == 1
      if check
         @assert is_subsymbol(toreplace.symbols[1], W.symbols[index])
      end

      first = change_pow(W.symbols[index],
         W.symbols[index].pow - toreplace.symbols[1].pow)
      last = change_pow(W.symbols[index], 0)

   else
      if check
         @assert is_subsymbol(toreplace.symbols[1], W.symbols[index])
         @assert W.symbols[index+1:index+n-2] == toreplace.symbols[2:end-1]
         @assert is_subsymbol(toreplace.symbols[end], W.symbols[index+n-1])
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

function replace_all!(W::GWord{T},subst_dict::Dict{GWord{T},GWord{T}}) where {T}
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

function replace_all(W::GWord{T},subst_dict::Dict{GWord{T},GWord{T}}) where {T}
   W = deepcopy(W)
   replace_all!(W, subst_dict)
   return W
end

###############################################################################
#
#   Misc
#
###############################################################################

function generate_balls(S::Vector{T}, Id::T=parent(first(S))(); radius=2, op=*) where T<:GWord
    sizes = Int[]
    B = [Id]
    for i in 1:radius
        BB = [op(i,j) for (i,j) in Base.product(B,S)]
        B = unique([B; vec(BB)])
        push!(sizes, length(B))
    end
    return B, sizes
end

function generate_balls{T<:RingElem}(S::Vector{T}, Id::T=one(parent(first(S))); radius=2, op=*)
    sizes = Int[]
    B = [Id]
    for i in 1:radius
        BB = [op(i,j) for (i,j) in Base.product(B,S)]
        B = unique([B; vec(BB)])
        push!(sizes, length(B))
    end
    return B, sizes
end

###############################################################################
#
#   Includes
#
###############################################################################

include("FreeGroup.jl")
include("FPGroups.jl")
include("AutGroup.jl")

include("DirectProducts.jl")
include("WreathProducts.jl")

end # of module Groups
