

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

end

FPGroup() = FPGroup(Vector{FPSymbol}(), Vector{GWord{FPSymbol}}())

function show(io::IO, G::FPGroup)
    print(io, "Finitely presented group on $(length(G.gens)) gens and $(length(G.rels)) relations:\n")
    print(io, "gens:\t", join([g.gen for g in G.gens], ","),"\n")
    print(io, "rels:\t", join([rel for rel in G.rels], ","),"\n")
end

function add_gen!{T<:GSymbol}(G::FPGroup, g::T)
    if !(g in G.gens)
        push!(G.gens, g)
    end
    return G
end

function add_rel!{T<:FPSymbol}(G::FPGroup, w::GWord{T})
    if !(w in G.rels)
        push!(G.rels, w)
    end
    return G
end

end #of module FinitelyPresentedGroups
