###############################################################################
#
#   FreeSymbol/FreeGroupElem/FreeGroup definition
#
###############################################################################

struct FreeSymbol <: GSymbol
   id::Symbol
   pow::Int
end

FreeGroupElem = GroupWord{FreeSymbol}

mutable struct FreeGroup <: AbstractFPGroup
   gens::Vector{FreeSymbol}

   function FreeGroup(gens::Vector{T}) where {T<:GSymbol}
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
###############################################################################

elem_type(::Type{FreeGroup}) = FreeGroupElem

parent_type(::Type{FreeGroupElem}) = FreeGroup

###############################################################################
#
#   FreeSymbol constructors
#
###############################################################################

FreeSymbol(s::Symbol) = FreeSymbol(s,1)
FreeSymbol(s::String) = FreeSymbol(Symbol(s))

FreeGroup(n::Int, symbol::String="f") = FreeGroup([Symbol(symbol,i) for i in 1:n])

FreeGroup(a::Vector) = FreeGroup(FreeSymbol.(a))

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (G::FreeGroup)(w::GroupWord{FreeSymbol})
   if length(syllables(w)) > 0
      for s in w.symbols
         i = findfirst(g -> g.id == s.id, G.gens)
         i == 0 && throw(DomainError(
            "Symbol $s does not belong to $G."))
         s.pow % G.gens[i].pow == 0 || throw(DomainError(
            "Symbol $s doesn't belong to $G."))
      end
   end
   w.parent = G
   return w
end

(G::FreeGroup)(s::FreeSymbol) = G(FreeGroupElem(s))

###############################################################################
#
#   Basic manipulation
#
###############################################################################

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::FreeGroup)
   print(io, "Free group on $(length(G.gens)) generators: ")
   join(io, G.gens, ", ")
end

###############################################################################
#
#   Comparison
#
###############################################################################

###############################################################################
#
#   Inversion
#
###############################################################################
