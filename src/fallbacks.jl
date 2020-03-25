# workarounds
Base.one(G::Generic.PermGroup) = Generic.Perm(G.n)
Base.one(r::NCRingElem) = one(parent(r))

# fallback definitions
# note: the user should implement those on type, when possible
Base.eltype(w::GW) where GW<:GWord = eltype(GW)
AbstractAlgebra.elem_type(G::Gr) where Gr <:AbstractFPGroup = elem_type(Gr)

AbstractAlgebra.parent_type(g::Gw) where Gw <:GWord = parent_type(parent(Gr))

function Base.one(G::Gr) where Gr <: AbstractFPGroup
    El = elem_type(G)
    id = El(eltype(El)[])
    id.parent = G
    return id
end
