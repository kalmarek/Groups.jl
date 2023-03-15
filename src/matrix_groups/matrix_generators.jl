struct MatrixElt{N,T,N²} <: Groups.GSymbol
    id::Symbol
    inv::Bool
    mat::StaticArrays.SMatrix{N,N,T,N²}

    function MatrixElt{N,T}(
        id::Symbol,
        mat::AbstractMatrix,
        inv::Bool = false,
    ) where {N,T}
        n = LinearAlgebra.checksquare(mat)
        @assert N == n
        @assert !iszero(LinearAlgebra.det(mat))
        return new{N,T,N^2}(id, inv, mat)
    end
end

function MatrixElt{N}(
    id::Symbol,
    mat::AbstractMatrix,
    inv::Bool = false,
) where {N}
    return MatrixElt{N,eltype(mat)}(id, mat, inv)
end

Base.show(io::IO, m::MatrixElt) = print(io, m.id, m.inv ? "⁻¹" : "")

Base.:(==)(m::MatrixElt, n::MatrixElt) = m.mat == n.mat

Base.hash(m::MatrixElt, h::UInt) = hash(m.mat, hash(typeof(m), h))

function Base.inv(m::MatrixElt{N,T}) where {N,T}
    return MatrixElt{N,T}(m.id, round.(T, inv(m.mat)), !m.inv)
end

matrix(m::MatrixElt) = m.mat
