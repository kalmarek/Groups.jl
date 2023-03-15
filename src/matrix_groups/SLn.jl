include("eltary_matrices.jl")

struct SpecialLinearGroup{N,T,R,S} <: AbstractMatrixGroup{N,T}
    base_ring::R
    alphabet::Alphabet{S}
    gens::Vector{S}

    function SpecialLinearGroup{N}(base_ring) where {N}
        S = [
            ElementaryMatrix{N}(i, j, one(base_ring)) for i in 1:N for
            j in 1:N if i ≠ j
        ]
        alphabet = Alphabet(S)

        T = eltype(base_ring)
        R = typeof(base_ring)
        St = eltype(S)

        return new{N,T,R,St}(base_ring, alphabet, S)
    end
end

GroupsCore.ngens(SL::SpecialLinearGroup) = length(SL.gens)

function Base.show(io::IO, ::SpecialLinearGroup{N,T}) where {N,T}
    return print(io, "SL{$N,$T}")
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    SL::SpecialLinearGroup{N,T},
) where {N,T}
    return print(io, "special linear group of $N×$N matrices over $T")
end
