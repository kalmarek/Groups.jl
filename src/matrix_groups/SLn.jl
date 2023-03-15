include("eltary_matrices.jl")

struct SpecialLinearGroup{N,T,R,A,S} <: AbstractMatrixGroup{N,T}
    base_ring::R
    alphabet::A
    gens::S

    function SpecialLinearGroup{N}(base_ring) where {N}
        S = [ElementaryMatrix{N}(i, j, one(base_ring)) for i in 1:N for j in 1:N if i ≠ j]
        alphabet = Alphabet(S)

        return new{
            N,
            eltype(base_ring),
            typeof(base_ring),
            typeof(alphabet),
            typeof(S)
        }(base_ring, alphabet, S)
    end
end

GroupsCore.ngens(SL::SpecialLinearGroup{N}) where {N} = N^2 - N

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
