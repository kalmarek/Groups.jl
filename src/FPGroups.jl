###############################################################################
#
#   FPSymbol/FPGroupElem/FPGroup definition
#
###############################################################################

immutable FPSymbol <: GSymbol
   str::String
   pow::Int
end

typealias FPGroupElem GWord{FPSymbol}

type FPGroup <: Group
   gens::Vector{FPSymbol}
   #     order::Vector{T}
   #     fastmult_table::Array{Int,2}
   function FPGroup{T<:GSymbol}(gens::Vector{T})
      G = new(gens)
      G.gens = gens
      return G
   end
end

export FPSymbol, FPGroupElem, FPGroup, generators

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

parent_type(::Type{FPGroupElem}) = FPGroup

elem_type(::FPGroup) = FPGroupElem

###############################################################################
#
#   FPSymbol constructors
#
###############################################################################

FPSymbol(s::String) = FPSymbol(s,1)

FPGroup(a::Vector{String}) = FPGroup([FPSymbol(i) for i in a])

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (G::FPGroup)()
   id = FPGroupElem(FPSymbol("", 0))
   id.parent = G
   return id
end

function (G::FPGroup)(w::GWord)
   eltype(w.symbols) == FPSymbol || throw("Can not coerce $w to FPGroup $G.")
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

(G::FPGroup)(s::FPSymbol) = G(FPGroupElem(s))

###############################################################################
#
#   Basic manipulation
#
###############################################################################

hash(s::FPSymbol, h::UInt) = hash(s.str, hash(s.pow, hash(FPSymbol, h)))

change_pow(s::FPSymbol, n::Int) = FPSymbol(s.str, n)

length(s::FPSymbol) = abs(s.pow)

generators(G::FPGroup) = [G(FPGroupElem(g)) for g in G.gens]

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::FPGroup)
   print(io, "Finitely presented group on $(length(G.gens)) generators:\n")
   print(io, "gens:\t", join([g.str for g in G.gens], ", "))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function (==)(s::FPSymbol, t::FPSymbol)
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

inv(s::FPSymbol) = change_pow(s, -s.pow)

###############################################################################
#
#   Misc
#
###############################################################################

# function add_rel!{T<:FPSymbol}(G::FPGroup, w::GWord{T})
#    if !(w in G.rels)
#       w = G(w)
#       push!(G.rels, w)
#    end
#    return G
# end
#
# function quotientgroup(G::FPGroup, rels::Vector{FPGroupElem})
#    for r in rels
#       parent(r) == G || throw("Can not form quotient group: $r is not an element of $G")
#    end
#    H = deepcopy(G)
#    for rel in rels
#       add_rel!(H, rel)
#    end
#    return H
# end
