###############################################################################
#
#   FreeSymbol/FreeGroupElem/FreeGroup definition
#
###############################################################################

immutable FreeSymbol <: GSymbol
   str::String
   pow::Int
end

typealias FreeGroupElem GWord{FreeSymbol}

type FreeGroup <: FPGroup
   gens::Vector{FreeSymbol}
   #     order::Vector{T}
   #     fastmult_table::Array{Int,2}
   function FreeGroup{T<:GSymbol}(gens::Vector{T})
      G = new(gens)
      G.gens = gens
      return G
   end
end

export FreeGroupElem, FreeGroup, generators

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent_type(::Type{FreeGroupElem}) = FreeGroup

elem_type(::FreeGroup) = FreeGroupElem

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
   id = FreeGroupElem(FreeSymbol("", 0))
   id.parent = G
   return id
end

function (G::FreeGroup)(w::GWord)
   eltype(w.symbols) == FreeSymbol || throw("Can not coerce $w to FreeGroup $G.")
   if length(w) > 0
      for s in w.symbols
         i = findfirst(g -> g.str == s.str, G.gens)
         i == 0 && throw("Symbol $s does not belong to $G.")
         s.pow % G.gens[i].pow == 0 || throw("Symbol $s doesn't belong to $G.")
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

generators(G::FreeGroup) = [G(FreeGroupElem(g)) for g in G.gens]

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::FreeGroup)
   print(io, "Finitely presented group on $(length(G.gens)) generators:\n")
   print(io, "gens:\t", join([g.str for g in G.gens], ", "))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function (==)(s::FreeSymbol, t::FreeSymbol)
   isone(s) && isone(t) && return true
   s.str == t.str || return false
   s.pow == t.pow || return false
   return true
end

###############################################################################
#
#   Inversion
#
###############################################################################

inv(s::FreeSymbol) = change_pow(s, -s.pow)

###############################################################################
#
#   Misc
#
###############################################################################

# function add_rel!{T<:FreeSymbol}(G::FreeGroup, w::GWord{T})
#    if !(w in G.rels)
#       w = G(w)
#       push!(G.rels, w)
#    end
#    return G
# end
#
# function quotientgroup(G::FreeGroup, rels::Vector{FreeGroupElem})
#    for r in rels
#       parent(r) == G || throw("Can not form quotient group: $r is not an element of $G")
#    end
#    H = deepcopy(G)
#    for rel in rels
#       add_rel!(H, rel)
#    end
#    return H
# end
