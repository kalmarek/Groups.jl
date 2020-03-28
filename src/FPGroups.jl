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
    rels::Dict{FreeGroupElem, FreeGroupElem}

    function FPGroup(gens::Vector{T}, rels::Dict{FreeGroupElem, FreeGroupElem}) where {T<:GSymbol}
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

AbstractAlgebra.elem_type(::Type{FPGroup}) = FPGroupElem
AbstractAlgebra.parent_type(::Type{FPGroupElem}) = FPGroup

###############################################################################
#
#   FPSymbol constructors
#

FPSymbol(s::Symbol) = FPSymbol(s, 1)
FPSymbol(s::String) = FPSymbol(Symbol(s))
FPSymbol(s::GSymbol) = FPSymbol(s.id, s.pow)

FPGroup(n::Int, symbol::String="f") = FPGroup([Symbol(symbol,i) for i in 1:n])
FPGroup(a::AbstractVector) = FPGroup([FPSymbol(i) for i in a])
FPGroup(gens::Vector{FPSymbol}) = FPGroup(gens, Dict{FreeGroupElem, FreeGroupElem}())

FPGroup(H::FreeGroup) = FPGroup([FPSymbol(s) for s in H.gens])

###############################################################################
#
#   Parent object call overloads
#

function (G::FPGroup)(w::GWord)
    if isempty(w)
        return one(G)
    end

    @boundscheck for s in syllables(w)
        i = findfirst(g -> g.id == s.id, G.gens)
        i == 0 && throw(DomainError("Symbol $s does not belong to $G."))
        s.pow % G.gens[i].pow != 0 && throw(
            DomainError("Symbol $s doesn't belong to $G."))
    end

    w = FPGroupElem(FPSymbol.(syllables(w)))
    setparent!(w, G)
    return reduce!(w)
end

(G::FPGroup)(s::GSymbol) = G(FPGroupElem(s))

###############################################################################
#
#   String I/O
#

function show(io::IO, G::FPGroup)
    print(io, "FPgroup on $(length(G.gens)) generators ")
    strrels = join(G.rels, ", ")
    if length(strrels) > 200
       print(io, "⟨ ", join(G.gens, ", "), " | $(length(G.rels)) relation(s) ⟩.")
    else
       print(io, "⟨ ", join(G.gens, ", "), " | ", join(G.rels, ", "), " ⟩.")
    end
end

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

freepreimage(G::FPGroup) = parent(first(keys(G.rels)))
freepreimage(g::FPGroupElem) = freepreimage(parent(g))(syllables(g))

function add_rels!(G::FPGroup, newrels::Dict{FreeGroupElem,FreeGroupElem})
    for w in keys(newrels)
        haskey(G.rels, w) && continue
        G.rels[w] = newrels[w]
    end
    return G
end

function Base.:/(G::FPGroup, newrels::Vector{FPGroupElem})
    for r in newrels
        parent(r) == G || throw(DomainError(
        "Can not form quotient group: $r is not an element of $G"))
    end
    H = deepcopy(G)
    F = freepreimage(H)
    newrels = Dict(freepreimage(r) => one(F) for r in newrels)
    add_rels!(H, newrels)
    return H
end

function Base.:/(F::FreeGroup, rels::Vector{FreeGroupElem})
    for r in rels
        parent(r) == F || throw(DomainError(
         "Can not form quotient group: $r is not an element of $F"))
    end
    G = FPGroup(FPSymbol.(F.gens))
    G.rels = Dict(rel => one(F) for rel in unique(rels))
    return G
end
