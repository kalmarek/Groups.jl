struct Transvection <: GSymbol
    id::Symbol
    ij::UInt8
    inv::Bool

    function Transvection(id::Symbol, i::Integer, j::Integer, inv = false)
        @assert id in (:ϱ, :λ)
        return new(id, _indices(UInt8(i), UInt8(j)), inv)
    end
end

ϱ(i, j) = Transvection(:ϱ, i, j)
λ(i, j) = Transvection(:λ, i, j)

_indices(ij::UInt8) = (ij & 0xf0) >> 4, (ij & 0x0f)

function _indices(i::UInt8, j::UInt8)
    @boundscheck @assert i < typemax(i) ÷ 2
    @boundscheck @assert j < typemax(j) ÷ 2
    sizeof
    return (i << 4) + j
end

indices(t::Transvection) = Int.(_indices(t.ij))

function Base.getproperty(t::Transvection, s::Symbol)
    s === :i && return first(indices(t))
    s === :j && return last(indices(t))
    return Core.getfield(t, s)
end

function Base.show(io::IO, t::Transvection)
    id = if t.id === :ϱ
        'ϱ'
    else # if t.id === :λ
        'λ'
    end
    print(io, id, subscriptify(t.i), '.', subscriptify(t.j))
    t.inv && print(io, "^-1")
end

Base.inv(t::Transvection) = Transvection(t.id, _indices(t.ij)..., !t.inv)

Base.:(==)(t::Transvection, s::Transvection) =
    t.id === s.id && t.ij == s.ij && t.inv == s.inv
Base.hash(t::Transvection, h::UInt) = hash(t.id, hash(t.ij, hash(t.inv, h)))

Base.@propagate_inbounds function evaluate!(v::NTuple{T, N}, t::Transvection, A::Alphabet, tmp=one(first(v))) where {T, N}
    i, j = indices(t)
    @assert i ≤ length(v) && j ≤ length(v)

    @inbounds begin
        if t.id === :ϱ
            if !t.inv
                append!(word(v[i]), word(v[j]))
            else
                # append!(word(v[i]), inv(A, word(v[j])))
                for l in Iterators.reverse(word(v[j]))
                    push!(word(v[i]), inv(A, l))
                end
            end
        else # if t.id === :λ
            if !t.inv
                # prepend!(word(v[i]), word(v[j]))
                for l in Iterators.reverse(word(v[j]))
                    pushfirst!(word(v[i]), l)
                end
            else
                # prepend!(word(v[i]), inv(A, word(v[j])))
                for l in word(v[j])
                    pushfirst!(word(v[i]), inv(A, l))
                end
            end
        end

        _setnormalform!(v[i], false)
        _setvalidhash!(v[i], false)
    end
    normalform!(tmp, v[i])
    copyto!(v[i], tmp)

    return v
end

struct PermRightAut <: GSymbol
    perm::Vector{UInt8}

    function PermRightAut(p::AbstractVector{<:Integer})
        @assert sort(p) == 1:length(p)
        return new(p)
    end
end

function Base.show(io::IO, p::PermRightAut)
    print(io, 'σ')
    join(io, (subscriptify(Int(i)) for i in p.perm))
end

Base.inv(p::PermRightAut) = PermRightAut(invperm(p.perm))

Base.:(==)(p::PermRightAut, q::PermRightAut) = p.perm == q.perm
Base.hash(p::PermRightAut, h::UInt) = hash(p.perm, hash(PermRightAut, h))

function evaluate!(v::NTuple{T, N}, p::PermRightAut, ::Alphabet, tmp=one(first(v))) where {T, N}
    return v[p.perm]
end
