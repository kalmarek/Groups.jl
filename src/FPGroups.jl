

immutable FPSymbol <: GSymbol
   str::String
   pow::Int
end

typealias FPGroupElem GWord{FPSymbol}

type FPGroup <: Group
   gens::Vector{FPSymbol}
   rels::Vector{FPGroupElem}
   #     order::Vector{T}
   #     fastmult_table::Array{Int,2}
   function FPGroup{T<:GSymbol}(gens::Vector{T}, rels::Vector{GWord{T}})
      G = new(gens, rels)
      G.gens = gens
      rels = [G(r) for r in rels]
      G.rels = rels
      return G
   end
end

export FPSymbol, FPGroupElem, FPGroup, generators

parent_type(::Type{FPGroupElem}) = FPGroup

elem_type(::FPGroup) = FPGroupElem


FPSymbol(s::String) = FPSymbol(s,1)

FPGroup(a::Vector{String}) = FPGroup([FPSymbol(i) for i in a], FPGroupElem[])

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


hash(s::FPSymbol, h::UInt) = hash(s.str, hash(s.pow, hash(FPSymbol, h)))

isone(s::FPSymbol) = s.pow == 0

change_pow(s::FPSymbol, n::Int) = FPSymbol(s.str, n)

length(s::FPSymbol) = abs(s.pow)

generators(G::FPGroup) = [G(FPGroupElem(g)) for g in G.gens]

function show(io::IO, G::FPGroup)
    print(io, "Finitely presented group on $(length(G.gens)) gens and $(length(G.rels)) relations:\n")
    print(io, "gens:\t", join([g.gen for g in G.gens], ","),"\n")
    print(io, "rels:\t", join([rel for rel in G.rels], ","),"\n")
end

function show(io::IO, s::FPSymbol)
   if isone(s)
      print(io, "(id)")
   elseif s.pow == 1
      print(io, s.str)
   else
      print(io, (s.str)*"^$(s.pow)")
   end
end

function (==)(s::FPSymbol, t::FPSymbol)
   isone(s) && isone(t) && return true
   s.str == t.str || return false
   s.pow == t.pow || return false
   return true
end




inv(s::FPSymbol) = change_pow(s, -s.pow)

function add_rel!{T<:FPSymbol}(G::FPGroup, w::GWord{T})
    if !(w in G.rels)
        push!(G.rels, w)
    end
    return G
end

end #of module FinitelyPresentedGroups
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
