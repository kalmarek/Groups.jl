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
include("arithmetic.jl")
include("findreplace.jl")

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
    else
        G = parent(W)
        w = T([inv(s) for s in Iterators.reverse(syllables(W))])
        return G(w)
    end
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
