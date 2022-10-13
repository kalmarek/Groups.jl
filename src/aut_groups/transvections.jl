struct Transvection <: GSymbol
    id::Symbol
    i::UInt16
    j::UInt16
    inv::Bool

    function Transvection(id::Symbol, i::Integer, j::Integer, inv = false)
        @assert id in (:ϱ, :λ)
        return new(id, i, j, inv)
    end
end

ϱ(i, j) = Transvection(:ϱ, i, j)
λ(i, j) = Transvection(:λ, i, j)

function Base.show(io::IO, t::Transvection)
    id = if t.id === :ϱ
        'ϱ'
    else # if t.id === :λ
        'λ'
    end
    print(io, id, subscriptify(t.i), '.', subscriptify(t.j))
    t.inv && print(io, "^-1")
end

Base.inv(t::Transvection) = Transvection(t.id, t.i, t.j, !t.inv)

Base.:(==)(t::Transvection, s::Transvection) =
    t.id === s.id && t.i == s.i && t.j == s.j && t.inv == s.inv

Base.hash(t::Transvection, h::UInt) = hash(hash(t.id, hash(t.i)), hash(t.j, hash(t.inv, h)))

Base.@propagate_inbounds @inline function evaluate!(
    v::NTuple{T,N},
    t::Transvection,
    tmp = one(first(v)),
) where {T,N}
    i, j = t.i, t.j
    @assert 1 ≤ i ≤ length(v) && 1 ≤ j ≤ length(v)

    A = alphabet(parent(first(v)))

    @inbounds begin
        if t.id === :ϱ
            if !t.inv
                append!(word(v[i]), word(v[j]))
            else
                # append!(word(v[i]), inv(word(v[j]), A))
                for l in Iterators.reverse(word(v[j]))
                    push!(word(v[i]), inv(l, A))
                end
            end
        else # if t.id === :λ
            if !t.inv
                # prepend!(word(v[i]), word(v[j]))
                for l in Iterators.reverse(word(v[j]))
                    pushfirst!(word(v[i]), l)
                end
            else
                # prepend!(word(v[i]), inv(word(v[j]), A))
                for l in word(v[j])
                    pushfirst!(word(v[i]), inv(l, A))
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

evaluate!(v::NTuple{T,N}, p::PermRightAut, tmp = nothing) where {T,N} = v[p.perm]
