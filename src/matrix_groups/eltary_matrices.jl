struct ElementaryMatrix{N,T} <: Groups.GSymbol
    i::Int
    j::Int
    val::T
    ElementaryMatrix{N}(i, j, val=1) where {N} =
        (@assert i ≠ j; new{N,typeof(val)}(i, j, val))
end

function Base.show(io::IO, e::ElementaryMatrix)
    print(io, 'E', Groups.subscriptify(e.i), Groups.subscriptify(e.j))
    !isone(e.val) && print(io, "^$(e.val)")
end

Base.:(==)(e::ElementaryMatrix{N}, f::ElementaryMatrix{N}) where {N} =
    e.i == f.i && e.j == f.j && e.val == f.val

Base.hash(e::ElementaryMatrix, h::UInt) =
    hash(typeof(e), hash((e.i, e.j, e.val), h))

Base.inv(e::ElementaryMatrix{N}) where {N} =
    ElementaryMatrix{N}(e.i, e.j, -e.val)

function matrix(e::ElementaryMatrix{N,T}) where {N,T}
    m = StaticArrays.MMatrix{N,N,T}(LinearAlgebra.I)
    m[e.i, e.j] = e.val
    x = StaticArrays.SMatrix{N,N}(m)
    return x
end
