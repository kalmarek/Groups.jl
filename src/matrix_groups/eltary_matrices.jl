struct ElementaryMatrix{N,T} <: Groups.GSymbol
    i::Int
    j::Int
    val::T
    function ElementaryMatrix{N}(i, j, val = 1) where {N}
        return (@assert i ≠ j; new{N,typeof(val)}(i, j, val))
    end
end

function Base.show(io::IO, e::ElementaryMatrix)
    print(io, 'E', Groups.subscriptify(e.i), Groups.subscriptify(e.j))
    return !isone(e.val) && print(io, "^$(e.val)")
end

function Base.:(==)(e::ElementaryMatrix{N}, f::ElementaryMatrix{N}) where {N}
    return e.i == f.i && e.j == f.j && e.val == f.val
end

function Base.hash(e::ElementaryMatrix, h::UInt)
    return hash(typeof(e), hash((e.i, e.j, e.val), h))
end

function Base.inv(e::ElementaryMatrix{N}) where {N}
    return ElementaryMatrix{N}(e.i, e.j, -e.val)
end

function matrix(e::ElementaryMatrix{N,T}) where {N,T}
    m = StaticArrays.MMatrix{N,N,T}(LinearAlgebra.I)
    m[e.i, e.j] = e.val
    x = StaticArrays.SMatrix{N,N}(m)
    return x
end
