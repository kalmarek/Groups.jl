###############################################################################
#
#   FreeSymbol/FreeGroupElem/FreeGroup definition
#

struct FreeSymbol <: GSymbol
   id::Symbol
   pow::Int
end

FreeGroupElem = GroupWord{FreeSymbol}

mutable struct FreeGroup <: AbstractFPGroup
   gens::Vector{FreeSymbol}

   function FreeGroup(gens::AbstractVector{T}) where {T<:GSymbol}
      G = new(gens)
      G.gens = gens
      return G
   end
end

export FreeGroupElem, FreeGroup

###############################################################################
#
#   Type and parent object methods
#

Base.eltype(::Type{FreeGroup}) = FreeGroupElem
GroupsCore.parent_type(::Type{FreeGroupElem}) = FreeGroup

###############################################################################
#
#   FreeSymbol constructors
#

FreeSymbol(s::Symbol) = FreeSymbol(s,1)
FreeSymbol(s::AbstractString) = FreeSymbol(Symbol(s))
FreeSymbol(s::GSymbol) = FreeSymbol(s.id, s.pow)

FreeGroup(n::Int, symbol::String="f") = FreeGroup([Symbol(symbol,i) for i in 1:n])
FreeGroup(a::AbstractVector) = FreeGroup(FreeSymbol.(a))

###############################################################################
#
#   Parent object call overloads
#

function (G::FreeGroup)(w::GroupWord{FreeSymbol})
   for s in syllables(w)
      i = findfirst(g -> g.id == s.id, G.gens)
      isnothing(i) && throw(DomainError(
         "Symbol $s does not belong to $G."))
      s.pow % G.gens[i].pow == 0 || throw(DomainError(
         "Symbol $s doesn't belong to $G."))
   end
   setparent!(w, G)
   return reduce!(w)
end

(G::FreeGroup)(s::GSymbol) = G(FreeGroupElem(s))
(G::FreeGroup)(v::AbstractVector{<:GSymbol}) = G(FreeGroupElem(FreeSymbol.(v)))

###############################################################################
#
#   String I/O
#

function Base.show(io::IO, G::FreeGroup)
   print(io, "Free group on $(length(G.gens)) generators: ")
   join(io, G.gens, ", ")
end
