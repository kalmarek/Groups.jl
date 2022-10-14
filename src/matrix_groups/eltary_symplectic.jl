struct ElementarySymplectic{N,T} <: Groups.GSymbol
    symbol::Symbol
    i::Int
    j::Int
    val::T
    function ElementarySymplectic{N}(s::Symbol, i::Integer, j::Integer, val=1) where {N}
        @assert s ∈ (:A, :B)
        @assert iseven(N)
        n = N ÷ 2
        if s === :A
            @assert 1 ≤ i ≤ n && 1 ≤ j ≤ n && i ≠ j
        elseif s === :B
            @assert xor(1 ≤ i ≤ n, 1 ≤ j ≤ n) && xor(n < i ≤ N, n < j ≤ N)
        end
        return new{N,typeof(val)}(s, i, j, val)
    end
end

function Base.show(io::IO, s::ElementarySymplectic)
    i, j = Groups.subscriptify(s.i), Groups.subscriptify(s.j)
    print(io, s.symbol, i, j)
    !isone(s.val) && print(io, "^$(s.val)")
end

_ind(s::ElementarySymplectic{N}) where {N} = (s.i, s.j)
_local_ind(N_half::Integer, i::Integer) = ifelse(i <= N_half, i, i - N_half)
function _dual_ind(s::ElementarySymplectic{N}) where {N}
    if s.symbol === :A && return _ind(s)
    else#if s.symbol === :B
        return _dual_ind(N ÷ 2, s.i, s.j)
    end
end

function _dual_ind(N_half, i, j)
    @assert i <= N_half < j || j <= N_half < i
    if i <= N_half # && j > N_half
        i, j = j - N_half, i + N_half
    else
        i, j = j + N_half, i - N_half
    end
    return i, j
end

function Base.:(==)(s::ElementarySymplectic{N}, t::ElementarySymplectic{M}) where {N,M}
    N == M || return false
    s.symbol == t.symbol || return false
    s.val == t.val || return false
    return _ind(t) == _ind(s) || _ind(t) == _dual_ind(s)
end

Base.hash(s::ElementarySymplectic, h::UInt) =
    hash(Set([_ind(s); _dual_ind(s)]), hash(s.symbol, hash(s.val, h)))

LinearAlgebra.transpose(s::ElementarySymplectic{N}) where {N} =
    ElementarySymplectic{N}(s.symbol, s.j, s.i, s.val)

Base.inv(s::ElementarySymplectic{N}) where {N} =
    ElementarySymplectic{N}(s.symbol, s.i, s.j, -s.val)

function matrix_repr(s::ElementarySymplectic{N,T}) where {N,T}
    @assert iseven(N)
    n = div(N, 2)
    m = StaticArrays.MMatrix{N,N,T}(LinearAlgebra.I)
    i, j = _ind(s)
    m[i, j] = s.val
    if s.symbol === :A
        m[n+j, n+i] = -s.val
    else#if s.symbol === :B
        if i > n
            m[j+n, i-n] = s.val
        else
            m[j-n, i+n] = s.val
        end
    end
    return StaticArrays.SMatrix{N,N}(m)
end
