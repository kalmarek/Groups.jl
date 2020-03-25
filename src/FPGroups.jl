###############################################################################
#
#   FPSymbol/FPGroupElem/FPGroup definition
#
###############################################################################

struct FPSymbol <: GSymbol
   id::Symbol
   pow::Int
end

FPGroupElem = GroupWord{FPSymbol}

mutable struct FPGroup <: AbstractFPGroup
   gens::Vector{FPSymbol}
   rels::Dict{FPGroupElem, FPGroupElem}

   function FPGroup(gens::Vector{T}, rels::Dict{FPGroupElem, FPGroupElem}) where {T<:GSymbol}
      G = new(gens)
      G.rels = Dict(G(k) => G(v) for (k,v) in rels)
      return G
   end
end

export FPGroupElem, FPGroup

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

elem_type(::Type{FPGroup}) = FPGroupElem
parent_type(::Type{FPGroupElem}) = FPGroup

###############################################################################
#
#   FPSymbol constructors
#
###############################################################################

FPSymbol(s::Symbol) = FPSymbol(s, 1)
FPSymbol(s::String) = FPSymbol(Symbol(s))
FPSymbol(s::GSymbol) = FPSymbol(s.id, s.pow)

convert(::Type{FPSymbol}, s::FreeSymbol) = FPSymbol(s.id, s.pow)

FPGroup(gens::Vector{FPSymbol}) = FPGroup(gens, Dict{FPGroupElem, FPGroupElem}())

FPGroup(a::Vector{String}) = FPGroup([FPSymbol(i) for i in a])

FPGroup(n::Int, symbol::String="f") = FPGroup(["$symbol$i" for i in 1:n])
FPGroup(H::FreeGroup) = FPGroup([FPSymbol(s) for s in H.gens])

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (G::FPGroup)(w::GWord)
   if isempty(w)
      return one(G)
   end

   if eltype(w.symbols) == FreeSymbol
      w = FPGroupElem(FPSymbol.(w.symbols))
   end

   if eltype(w.symbols) == FPSymbol
      for s in w.symbols
         i = findfirst(g -> g.id == s.id, G.gens)
         i == 0 && throw(DomainError(
            "Symbol $s does not belong to $G."))
         s.pow % G.gens[i].pow == 0 || throw(DomainError(
         "Symbol $s doesn't belong to $G."))
      end
   end
   w.parent = G
   return reduce!(w)
end

(G::FPGroup)(s::FPSymbol) = G(FPGroupElem(s))

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

function show(io::IO, G::FPGroup)
   print(io, "FPgroup on $(length(G.gens)) generators ")
   strrels = join(G.rels, ", ")
   if length(strrels) > 300
      print(io, "⟨ ", join(G.gens, ", "), " | $(length(G.rels)) relation(s) ⟩.")
   else
      print(io, "⟨ ", join(G.gens, ", "), " | ", join(G.rels, ", "), " ⟩.")
   end
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

###############################################################################
#
#   Binary operations
#
###############################################################################

function reduce!(W::FPGroupElem)
    reduced = false
    while !reduced
        W = replace(W, parent(W).rels)
        reduced = freereduce!(Bool, W)
    end
    return W
end

###############################################################################
#
#   Misc
#
###############################################################################

function add_rels!(G::FPGroup, newrels::Dict{FPGroupElem,FPGroupElem})
   for w in keys(newrels)
      if !(w in keys(G.rels))
         G.rels[w] = G(newrels[w])
      end
   end
end

function Base.:/(G::FPGroup, newrels::Vector{FPGroupElem})
   for r in newrels
      parent(r) == G || throw(DomainError(
      "Can not form quotient group: $r is not an element of $G"))
   end
   H = deepcopy(G)
   newrels = Dict(H(r) => one(H) for r in newrels)
   add_rels!(H, newrels)
   return H
end

function Base.:/(G::FreeGroup, rels::Vector{FreeGroupElem})
   for r in rels
      parent(r) == G || throw(DomainError(
         "Can not form quotient group: $r is not an element of $G"))
   end
   H = FPGroup(deepcopy(G))
   H.rels = Dict(H(rel) => one(H) for rel in unique(rels))
   return H
end
