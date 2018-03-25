###############################################################################
#
#   FPSymbol/FPGroupElem/FPGroup definition
#
###############################################################################

struct FPSymbol <: GSymbol
   str::String
   pow::Int
end

FPGroupElem = GroupWord{FPSymbol}

mutable struct FPGroup <: AbstractFPGroup
   gens::Vector{FPSymbol}
   rels::Dict{FPGroupElem, FPGroupElem}

   function FPGroup{T<:GSymbol}(gens::Vector{T}, rels::Dict{FPGroupElem, FPGroupElem})
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

parent_type(::Type{FPGroupElem}) = FPGroup

elem_type(::FPGroup) = FPGroupElem

###############################################################################
#
#   FPSymbol constructors
#
###############################################################################

FPSymbol(s::String) = FPSymbol(s,1)

convert(::Type{FPSymbol}, s::FreeSymbol) = FPSymbol(s.str, s.pow)

FPGroup(gens::Vector{FPSymbol}) = FPGroup(gens, Dict{FPGroupElem, FPGroupElem}())

FPGroup(a::Vector{String}) = FPGroup([FPSymbol(i) for i in a])

FPGroup(n::Int, symbol::String="f") = FPGroup(["$symbol$i" for i in 1:n])
FPGroup(H::FreeGroup) = FPGroup([FPSymbol(s) for s in H.gens])

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
   if length(w) == 0
      return G()
   end

   if eltype(w.symbols) == FreeSymbol
      w = FPGroupElem(w.symbols)
   end

   if eltype(w.symbols) == FPSymbol
      for s in w.symbols
         i = findfirst(g -> g.str == s.str, G.gens)
         i == 0 && throw("Symbol $s does not belong to $G.")
         s.pow % G.gens[i].pow == 0 || throw("Symbol $s doesn't belong to $G.")
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

hash(s::FPSymbol, h::UInt) = hash(s.str, hash(s.pow, hash(FPSymbol, h)))

change_pow(s::FPSymbol, n::Int) = FPSymbol(s.str, n)

length(s::FPSymbol) = abs(s.pow)

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
#   Binary operations
#
###############################################################################

(*)(W::FPGroupElem, Z::FPGroupElem) = r_multiply(W, Z.symbols)
(*)(W::FPGroupElem, s::FPSymbol) = r_multiply(W, [s])
(*)(s::FPSymbol, W::FPGroupElem) = l_multiply(W, [s])

function reduce!(W::FPGroupElem)
    if length(W) < 2
        deleteat!(W.symbols, find(x -> x.pow == 0, W.symbols))
    else
        reduced = false
        while !reduced
            reduced = free_reduce!(W) || replace_all!(W, parent(W).rels)
        end
    end

    W.savedhash = hash(W.symbols, hash(typeof(W)))
    W.modified = false
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

function /(G::FPGroup, newrels::Vector{FPGroupElem})
   for r in rels
      parent(r) == G || throw("Can not form quotient group: $r is not an element of $G")
   end
   H = deepcopy(G)
   newrels = Dict(H(r) => H() for r in newrels)
   add_rels!(H, newrels)
   return H
end

function /(G::FreeGroup, rels::Vector{FreeGroupElem})
   for r in rels
      parent(r) == G || throw("Can not form quotient group: $r is not an element of $G")
   end
   H = FPGroup(G)
   H.rels = Dict(H(rel) => H() for rel in unique(rels))
   return H
end
