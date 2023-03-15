include("matrix_generators.jl")

struct MatrixGroup{N,T,R,S} <: AbstractMatrixGroup{N,T}
    base_ring::R
    alphabet::Alphabet{S}
    gens::Vector{S}
end

function MatrixGroup{N}(
    gens::AbstractVector{<:AbstractMatrix{T}},
    base_ring = T,
) where {N,T}
    S = map(enumerate(gens)) do (i, mat)
        id = Symbol('m', Groups.subscriptify(i))
        return MatrixElt{N}(id, mat)
    end
    alphabet = Alphabet(S)

    R = typeof(base_ring)
    St = eltype(S)

    return MatrixGroup{N,T,R,St}(base_ring, alphabet, S)
end

GroupsCore.ngens(M::MatrixGroup) = length(M.gens)
