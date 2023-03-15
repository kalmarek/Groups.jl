abstract type AbstractMatrixGroup{N,T} <: Groups.AbstractFPGroup end
const MatrixGroupElement{N,T} =
    Groups.AbstractFPGroupElement{<:AbstractMatrixGroup{N,T}}

function Base.isone(g::MatrixGroupElement{N,T}) where {N,T}
    return isone(word(g)) || isone(matrix(g))
end

function Base.:(==)(
    m1::M1,
    m2::M2,
) where {M1<:MatrixGroupElement,M2<:MatrixGroupElement}
    parent(m1) === parent(m2) || return false
    word(m1) == word(m2) && return true
    return matrix(m1) == matrix(m2)
end

Base.size(::MatrixGroupElement{N}) where {N} = (N, N)
Base.size(::MatrixGroupElement{N}, d) where {N} = ifelse(d::Integer <= 2, N, 1)
Base.eltype(::MatrixGroupElement{N,T}) where {N,T} = T

# three structural assumptions about matrix groups
Groups.word(m::MatrixGroupElement) = m.word
Base.parent(m::MatrixGroupElement) = m.parent
Groups.alphabet(M::AbstractMatrixGroup) = M.alphabet
Groups.rewriting(M::AbstractMatrixGroup) = alphabet(M)

Base.hash(m::MatrixGroupElement, h::UInt) = hash(matrix(m), hash(parent(m), h))

function matrix(m::MatrixGroupElement{N,T}) where {N,T}
    if isone(word(m))
        return StaticArrays.SMatrix{N,N,T}(LinearAlgebra.I)
    end
    A = alphabet(parent(m))
    return prod(matrix(A[l]) for l in word(m))
end

function Base.convert(
    ::Type{M},
    m::MatrixGroupElement,
) where {M<:AbstractMatrix}
    return convert(M, matrix(m))
end
(M::Type{<:AbstractMatrix})(m::MatrixGroupElement) = convert(M, m)

function Base.rand(
    rng::Random.AbstractRNG,
    rs::Random.SamplerTrivial{<:AbstractMatrixGroup},
)
    Mgroup = rs[]
    S = gens(Mgroup)
    return prod(
        g -> rand(rng, Bool) ? g : inv(g),
        rand(rng, S, rand(rng, 1:30)),
    )
end

function Base.show(io::IO, M::AbstractMatrixGroup)
    g = gens(M, 1)
    N = size(g, 1)
    return print(io, "H ⩽ GL{$N,$(eltype(g))}")
end

function Base.show(io::IO, ::MIME"text/plain", M::AbstractMatrixGroup)
    N = size(gens(M, 1), 1)
    ng = GroupsCore.ngens(M)
    return print(
        io,
        "subgroup of $N×$N invertible matrices with $(ng) generators",
    )
end

function Base.show(
    io::IO,
    mat::Groups.AbstractFPGroupElement{<:AbstractMatrixGroup},
)
    return KnuthBendix.print_repr(io, word(mat), alphabet(mat))
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    mat::Groups.AbstractFPGroupElement{<:AbstractMatrixGroup{N}},
) where {N}
    Groups.normalform!(mat)
    KnuthBendix.print_repr(io, word(mat), alphabet(mat))
    println(io, " ∈ ", parent(mat))
    return Base.print_array(io, matrix(mat))
end
