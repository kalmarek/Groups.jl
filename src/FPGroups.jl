

immutable FPSymbol <: GSymbol
   str::String
   pow::Int
end

typealias FPGroupElem GWord{FPSymbol}

end

# FPSymbol(x::String, G::Group) = FPSymbol(x,1,G)
# FPSymbol(s::GSymbol, G::Group) = FPSymbol(s.gen, s.pow, G)


immutable FPGroup <: Group
    gens::Vector{FPSymbol}
    rels::Vector{FPGroupElem}
#     order::Vector{T}
#     fast_multable::Array{Int,2}
    function FPGroup{T<:GSymbol}(gens::Vector{T}, rels::GWord{T})
        G = new()
        gens = [FPSymbol{G}(g.gen, G) for g in gens]
        rels = []
    end
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
