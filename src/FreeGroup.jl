###############################################################################
#
#   FreeSymbol/FreeGroupElem/FreeGroup definition
#
###############################################################################

struct FreeSymbol <: GSymbol
   str::String
   pow::Int
end

FreeGroupElem = GroupWord{FreeSymbol}

mutable struct FreeGroup <: AbstractFPGroup
   gens::Vector{FreeSymbol}

   function FreeGroup{T<:GSymbol}(gens::Vector{T})
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

FreeSymbol(s::String) = FreeSymbol(s,1)

FreeGroup(n::Int, symbol::String="f") = FreeGroup(["$symbol$i" for i in 1:n])

FreeGroup(a::Vector{String}) = FreeGroup([FreeSymbol(i) for i in a])

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (G::FreeGroup)()
   id = FreeGroupElem(FreeSymbol[])
   id.parent = G
   return id
end

function (G::FreeGroup)(w::GroupWord{FreeSymbol})
   if length(w) > 0
      for s in w.symbols
         i = findfirst(g -> g.str == s.str, G.gens)
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

hash(s::FreeSymbol, h::UInt) = hash(s.str, hash(s.pow, hash(FreeSymbol, h)))

change_pow(s::FreeSymbol, n::Int) = FreeSymbol(s.str, n)

length(s::FreeSymbol) = abs(s.pow)

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

inv(s::FreeSymbol) = change_pow(s, -s.pow)
