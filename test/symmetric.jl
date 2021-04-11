import AbstractAlgebra
using GroupsCore

# disambiguation
GroupsCore.order(
    ::Type{I},
    G::AbstractAlgebra.Generic.SymmetricGroup,
) where {I<:Integer} = I(factorial(G.n))

# disambiguation
GroupsCore.order(
    ::Type{I},
    g::AbstractAlgebra.Generic.Perm,
) where {I<:Integer} =
    I(foldl(lcm, length(c) for c in AbstractAlgebra.cycles(g)))

# correct the AA length:
Base.length(G::AbstractAlgebra.Generic.SymmetricGroup) = order(Int, G)

# genuinely new methods:
Base.IteratorSize(::Type{<:AbstractAlgebra.AbstractPermutationGroup}) = Base.HasLength()

function GroupsCore.gens(G::AbstractAlgebra.Generic.SymmetricGroup{I}) where {I}
    a, b = one(G), one(G)
    circshift!(a.d, b.d, -1)
    b.d[1], b.d[2] = 2, 1
    return [a, b]
end

Base.deepcopy_internal(g::AbstractAlgebra.Generic.Perm, ::IdDict) =
    AbstractAlgebra.Generic.Perm(deepcopy(g.d), false)
