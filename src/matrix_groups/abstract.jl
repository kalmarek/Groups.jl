abstract type MatrixGroup{N, T} <: Groups.AbstractFPGroup end
const MatrixGroupElement{N, T} = Groups.AbstractFPGroupElement{<:MatrixGroup{N, T}}

Base.isone(g::MatrixGroupElement{N, T}) where {N, T} = matrix_repr(g) == I

function Base.:(==)(m1::M1, m2::M2) where {M1<:MatrixGroupElement, M2<:MatrixGroupElement}
    parent(m1) === parent(m2) || return false
    word(m1) == word(m2) && return true
    return matrix_repr(m1) == matrix_repr(m2)
end

Base.size(m::MatrixGroupElement{N}) where N = (N, N)
Base.eltype(m::MatrixGroupElement{N, T}) where {N, T} = T

# three structural assumptions about matrix groups
Groups.word(sl::MatrixGroupElement) = sl.word
Base.parent(sl::MatrixGroupElement) = sl.parent
Groups.alphabet(M::MatrixGroup) = M.alphabet
Groups.rewriting(M::MatrixGroup) = alphabet(M)

Base.hash(sl::MatrixGroupElement, h::UInt) =
    hash(matrix_repr(sl), hash(parent(sl), h))

Base.getindex(sl::MatrixGroupElement, i, j) = matrix_repr(sl)[i,j]
# Base.iterate(sl::MatrixGroupElement) = iterate(sl.elts)
# Base.iterate(sl::MatrixGroupElement, state) = iterate(sl.elts, state)

function matrix_repr(m::MatrixGroupElement{N, T}) where {N, T}
    isempty(word(m)) && return StaticArrays.SMatrix{N, N, T}(I)
    A = alphabet(parent(m))
    return prod(matrix_repr(A[l]) for l in word(m))
end

function Base.rand(
    rng::Random.AbstractRNG,
    rs::Random.SamplerTrivial{<:MatrixGroup},
    )
    Mgroup = rs[]
    S = gens(Mgroup)
    return prod(g -> rand(Bool) ? g : inv(g), rand(S, rand(1:30)))
end
