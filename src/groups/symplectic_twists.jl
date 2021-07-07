struct ΡΛ
    id::Symbol
    A::Alphabet
    N::Int
end

function Base.getindex(rl::ΡΛ, i::Integer, j::Integer)
    @assert 1 ≤ i ≤ rl.N
    @assert 1 ≤ j ≤ rl.N
    @assert i ≠ j
    @assert rl.id ∈ (:λ, :ϱ)
    rl.id == :λ && return Word([rl.A[λ(i, j)]])
    rl.id == :ϱ && return Word([rl.A[ϱ(i, j)]])
end

function Te_diagonal(λ::ΡΛ, ϱ::ΡΛ, i::Integer)
    @assert λ.N == ϱ.N
    @assert λ.id == :λ && ϱ.id == :ϱ

    N = λ.N
    @assert iseven(N)
    n = N ÷ 2
    j = i + 1
    @assert 1 <= i < n

    A = λ.A

    # comments are for i,j = 1,2
    g = one(Word(Int[]))
    g *= λ[n+j, n+i]                   # β ↦ α*β
    g *= λ[n+i, i] * inv(A, ϱ[n+i, j]) # α ↦ a*α*b^-1
    g *= inv(A, λ[n+j, n+i])           # β ↦ b*α^-1*a^-1*α*β
    g *= λ[j, n+i] * inv(A, λ[j, i])   # b ↦ α
    g *= inv(A, λ[j, n+i])             # b ↦ b*α^-1*a^-1*α
    g *= inv(A, ϱ[j, n+i]) * ϱ[j, i]   # b ↦ b*α^-1*a^-1*α*b*α^-1
    g *= ϱ[j, n+i]                     # b ↦ b*α^-1*a^-1*α*b*α^-1*a*α*b^-1
    return g
end

function Te_lantern(A::Alphabet, b₀::T, a₁::T, a₂::T, a₃::T, a₄::T, a₅::T) where {T}
    a₀ = (a₁ * a₂ * a₃)^4 * inv(A, b₀)
    X = a₄ * a₅ * a₃ * a₄
    b₁ = inv(A, X) * a₀ * X
    Y = a₂ * a₃ * a₁ * a₂
    return inv(A, Y) * b₁ * Y # b₂
end

Ta(λ::ΡΛ, i::Integer) = (@assert λ.id == :λ;
λ[λ.N÷2+i, i])
Tα(λ::ΡΛ, i::Integer) = (@assert λ.id == :λ;
inv(λ.A, λ[i, λ.N÷2+i]))

function Te(λ::ΡΛ, ϱ::ΡΛ, i, j)
    @assert i ≠ j
    i, j = i < j ? (i, j) : (j, i)

    @assert λ.N == ϱ.N
    @assert λ.A == ϱ.A
    @assert λ.id == :λ && ϱ.id == :ϱ

    @assert 1 ≤ i ≤ λ.N
    @assert 1 ≤ j ≤ λ.N

    if j == i + 1
        return Te_diagonal(λ, ϱ, i)
    else
        return Te_lantern(
            λ.A,
            Ta(λ, i + 1),
            Ta(λ, i),
            Tα(λ, i),
            Te(λ, ϱ, i, i + 1),
            Tα(λ, i + 1),
            Te(λ, ϱ, i + 1, j),
        )
    end
end

function mcg_twists(G::AutomorphismGroup{<:FreeGroup})

    @assert iseven(ngens(object(G)))
    genus = ngens(object(G)) ÷ 2

    genus < 3 && throw("Not Implemented: genus = $genus < 3")

    A = KnuthBendix.alphabet(G)

    λ = ΡΛ(:λ, A, 2genus)
    ϱ = ΡΛ(:ϱ, A, 2genus)

    Tas = [Ta(λ, i) for i in 1:genus]
    Tαs = [Tα(λ, i) for i in 1:genus]

    idcs = ((i, j) for i in 1:genus for j in i+1:genus)
    Tes = [Te(λ, ϱ, i, j) for (i, j) in idcs]

    return Tas, Tαs, Tes
end

struct SymplecticMappingClass{N,T} <: GSymbol
    id::Symbol # :A, :B
    i::UInt
    j::UInt
    minus::Bool
    inv::Bool
    images::NTuple{N,T}
    invimages::NTuple{N,T}
    gens_idcs::Dict{Int, Int}
end

function SymplecticMappingClass(
    Σ::SurfaceGroup,
    sautFn,
    id::Symbol,
    i::Integer,
    j::Integer;
    minus = false,
    inverse = false,
)
    @assert i > 0 && j > 0
    id === :A && @assert i ≠ j
    @assert 2genus(Σ) == ngens(object(sautFn))

    A = KnuthBendix.alphabet(sautFn)
    λ = ΡΛ(:λ, A, 2genus(Σ))
    ϱ = ΡΛ(:ϱ, A, 2genus(Σ))

    w = if id === :A
        Te(λ, ϱ, i, j) *
        inv(A, Ta(λ, i)) *
        Tα(λ, i) *
        Ta(λ, i) *
        inv(A, Te(λ, ϱ, i, j)) *
        inv(A, Tα(λ, i)) *
        inv(A, Ta(λ, j))
    elseif id === :B
        if !minus
            if i ≠ j
                x = Ta(λ, j) * inv(A, Ta(λ, i)) * Tα(λ, j) * Te(λ, ϱ, i, j)
                δ = x * Tα(λ, i) * inv(A, x)
                Tα(λ, i) * Tα(λ, j) * inv(A, δ)
            else
                inv(A, Tα(λ, i))
            end
        else
            if i ≠ j
                Ta(λ, i) * Ta(λ, j) * inv(A, Te(λ, ϱ, i, j))
            else
                Ta(λ, i)
            end
        end
    else
        throw("Type not recognized: $id")
    end

    g = sautFn(w)

    d = ntuple(i->gens(Σ, i), ngens(Σ))

    img = evaluate!(deepcopy(d), g)
    invim = evaluate!(d, inv(g))

    img, invim = inverse ? (invim, img) : (img, invim)

    gens_idcs = Dict(alphabet(Σ)[Σ.gens[i]] => i for i in 1:ngens(Σ))

    res = SymplecticMappingClass(id, UInt(i), UInt(j), minus, inverse, img, invim, gens_idcs)

    return res
end

function Base.show(io::IO, smc::SymplecticMappingClass)
    smc.minus && print(io, 'm')
    if  smc.i < 10 && smc.j < 10
        print(io, smc.id, subscriptify(smc.i), subscriptify(smc.j))
    else
        print(io, smc.id, subscriptify(smc.i), ".", subscriptify(smc.j))
    end
    smc.inv && print(io, "^-1")
end

function Base.inv(m::SymplecticMappingClass)
    return SymplecticMappingClass(m.id, m.i, m.j, m.minus, !m.inv, m.invimages, m.images, m.gens_idcs)
end

function evaluate!(
    t::NTuple{N,T},
    smc::SymplecticMappingClass,
    A::Alphabet,
    tmp = one(first(t)),
) where {N,T}
    img = smc.images

    for elt in t
        copyto!(tmp, elt)
        resize!(word(elt), 0)
        for idx in word(tmp)
            # @show idx
            k = if haskey(smc.gens_idcs, idx)
                img[smc.gens_idcs[idx]]
            else
                inv(img[smc.gens_idcs[inv(A, idx)]])
            end
            append!(word(elt), word(k))
        end
        _setnormalform!(elt, false)
        _setvalidhash!(elt, false)

        normalform!(tmp, elt)
        copyto!(elt, tmp)
    end

    return t
end
