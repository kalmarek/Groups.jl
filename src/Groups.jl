module Groups

using AbstractAlgebra
import AbstractAlgebra: Group, GroupElem, Ring
import AbstractAlgebra: parent, parent_type, elem_type
import AbstractAlgebra: order, gens, matrix_repr

import Base: length, ==, hash, show, convert, eltype, iterate
import Base: inv, reduce, *, ^, power_by_squaring
import Base: findfirst, findnext, replace
import Base: deepcopy_internal

using LinearAlgebra
using Markdown


include("types.jl")
include("gsymbols.jl")
include("fallbacks.jl")
include("words.jl")
include("hashing.jl")
include("freereduce.jl")

include("FreeGroup.jl")
include("FPGroups.jl")
include("AutGroup.jl")

include("DirectPower.jl")
include("WreathProducts.jl")


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

issubsymbol(s::GSymbol, t::GSymbol) =
    s.id == t.id && (0 ≤ s.pow ≤ t.pow || 0 ≥ s.pow ≥ t.pow)

function issubsymbol(s::FreeSymbol, w::GWord, sindex::Integer)
    @boundscheck 1 ≤ sindex ≤ syllablelength(w) || throw(BoundsError(w, sindex))
    return issubsymbol(s, syllables(w)[sindex])
end

function issubword(z::GWord, w::GWord, sindex::Integer)
    isempty(z) && return true
    @boundscheck 1 ≤ sindex ≤ syllablelength(w) || throw(BoundsError(w, sindex))
    n = syllablelength(z)
    n == 1 && return issubsymbol(first(syllables(z)), syllables(w)[sindex])

    lastindex = sindex + n - 1
    lastindex > syllablelength(w) && return false

    issubsymbol(first(z), syllables(w)[sindex]) || return false
    issubsymbol(syllables(z)[end], syllables(w)[lastindex]) || return false
    for (zidx, widx) in zip(2:n-1, sindex+1:lastindex-1)
        syllables(z)[zidx] == syllables(w)[widx] || return false
    end
    return true
end

"""doc
Find the first syllable index k>=i such that Z < syllables(W)[k:k+syllablelength(Z)-1]
"""
function findnext(subword::GWord, word::GWord, start::Integer)
    @boundscheck 1 ≤ start ≤ syllablelength(word) || throw(BoundsError(word, start))
    isempty(subword) && return start
    stop = syllablelength(word) - syllablelength(subword) +1

    for idx in start:1:stop
        issubword(subword, word, idx) && return idx
    end
    return nothing
end

function findnext(s::FreeSymbol, word::GWord, start::Integer)
    @boundscheck 1 ≤ start ≤ syllablelength(word) || throw(BoundsError(word, start))
    isone(s) && return start
    stop = syllablelength(word)

    for idx in start:1:stop
        issubsymbol(s, word, idx) && return idx
    end
    return nothing
end

function findprev(subword::GWord, word::GWord, start::Integer)
    @boundscheck 1 ≤ start ≤ syllablelength(word) || throw(BoundsError(word, start))
    isempty(subword) && return start
    stop = 1

    for idx in start:-1:1
        issubword(subword, word, idx) && return idx
    end
    return nothing
end

function findprev(s::FreeSymbol, word::GWord, start::Integer)
    @boundscheck 1 ≤ start ≤ syllablelength(word) || throw(BoundsError(word, start))
    isone(s) && return start
    stop = 1

    for idx in start:-1:stop
        issubsymbol(s, word, idx) && return idx
    end
    return nothing
end

findfirst(subword::GWord, word::GWord) = findnext(subword, word, 1)
findlast(subword::GWord, word::GWord) =
    findprev(subword, word, syllablelength(word)-syllablelength(subword)+1)

function replace!(out::GW, W::GW, lhs_rhs::Pair{GS, T}; count::Integer=typemax(Int)) where
        {GS<:GSymbol, T<:GWord, GW<:GWord}
    (count == 0 || isempty(W)) && return W
    count < 0 && throw(DomainError(count, "`count` must be non-negative."))

    lhs, rhs = lhs_rhs

    sW = syllables(W)
    sW_idx = 1
    r = something(findnext(lhs, W, sW_idx), 0)

    sout = syllables(out)
    resize!(sout, 0)
    sizehint!(sout, syllablelength(W))

    c = 0

    while !iszero(r)
        append!(sout, view(sW, sW_idx:r-1))
        a, b = divrem(sW[r].pow, lhs.pow)

        if b != 0
            push!(sout, change_pow(sW[r], b))
        end

        append!(sout, repeat(syllables(rhs), a))

        sW_idx = r+1
        sW_idx > syllablelength(W) && break

        r = something(findnext(lhs, W, sW_idx), 0)
        c += 1
        c == count && break
    end
    append!(sout, sW[sW_idx:end])
    return freereduce!(out)
end

function replace!(out::GW, W::GW, lhs_rhs::Pair{T, T}; count::Integer=typemax(Int)) where
    {GW<:GWord, T <: GWord}
    (count == 0 || isempty(W)) && return W
    count < 0 && throw(DomainError(count, "`count` must be non-negative."))

    lhs, rhs = lhs_rhs
    lhs_slen = syllablelength(lhs)
    lhs_slen == 1 && return replace!(out, W, first(syllables(lhs))=>rhs; count=count)

    sW = syllables(W)
    sW_idx = 1
    r = something(findnext(lhs, W, sW_idx), 0)

    sout = syllables(out)
    resize!(sout, 0)
    sizehint!(sout, syllablelength(W))

    c = 0

    while !iszero(r)
        append!(sout, view(sW, sW_idx:r-1))

        exp = sW[r].pow - first(syllables(lhs)).pow
        if exp != 0
            push!(sout, change_pow(sW[r], exp))
        end

        append!(sout, syllables(rhs))

        exp = sW[r+lhs_slen-1].pow - last(syllables(lhs)).pow
        if exp != 0
            push!(sout, change_pow(sW[r+lhs_slen-1], exp))
        end

        sW_idx = r+lhs_slen
        sW_idx > syllablelength(W) && break

        r = something(findnext(lhs, W, sW_idx), 0)
        c += 1
        c == count && break
    end

    # copy the rest
    append!(sout, sW[sW_idx:end])
    return freereduce!(out)
end

function replace(W::GW, lhs_rhs::Pair{T, T}; count::Integer=typemax(Int)) where
    {GW<:GWord, T <: GWord}
    return replace!(one(W), W, lhs_rhs; count=count)
end

function replace(W::GW, subst_dict::Dict{T,T}) where {GW<:GWord, T<:GWord}
    out = W
    for toreplace in reverse!(sort!(collect(keys(subst_dict)), by=length))
        replacement = subst_dict[toreplace]
        if length(toreplace) > length(out)
            continue
        end
        out = replace(out, toreplace=>replacement)
    end
    return out
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
