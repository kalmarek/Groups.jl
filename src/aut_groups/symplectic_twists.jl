struct ΡΛ
    id::Symbol
    A::Alphabet
    N::Int
end

function Base.getindex(rl::ΡΛ, i::Integer, j::Integer)
    @assert 1 ≤ i ≤ rl.N "Got $i > $(rl.N)"
    @assert 1 ≤ j ≤ rl.N "Got $j > $(rl.N)"
    @assert i ≠ j
    @assert rl.id ∈ (:λ, :ϱ)
    rl.id == :λ && return Word([rl.A[λ(i, j)]])
    rl.id == :ϱ && return Word([rl.A[ϱ(i, j)]])
end

function Te_diagonal(λ::Groups.ΡΛ, ϱ::Groups.ΡΛ, i::Integer)
    @assert λ.N == ϱ.N
    @assert λ.id == :λ && ϱ.id == :ϱ

    N = λ.N
    @assert iseven(N)
    A = λ.A
    n = N ÷ 2
    j = i + 1

    if i == n
        τ = rotation_element(λ, ϱ)
        return inv(τ, A) * Te_diagonal(λ, ϱ, 1) * τ
    end

    @assert 1 <= i < n

    NI = (2n - 2i) + 1
    NJ = (2n - 2j) + 1
    I = (2n - 2i) + 2
    J = (2n - 2j) + 2

    g = one(Word(Int[]))
    g *= λ[NJ, NI]                    # β ↦ α*β
    g *= λ[NI, I] * inv(ϱ[NI, J], A)  # α ↦ a*α*b^-1
    g *= inv(λ[NJ, NI], A)            # β ↦ b*α^-1*a^-1*α*β
    g *= λ[J, NI] * inv(λ[J, I], A)   # b ↦ α
    g *= inv(λ[J, NI], A)             # b ↦ b*α^-1*a^-1*α
    g *= inv(ϱ[J, NI], A) * ϱ[J, I]   # b ↦ b*α^-1*a^-1*α*b*α^-1
    g *= ϱ[J, NI]                     # b ↦ b*α^-1*a^-1*α*b*α^-1*a*α*b^-1

    return g
end

function Te_lantern(A::Alphabet, b₀::T, a₁::T, a₂::T, a₃::T, a₄::T, a₅::T) where {T}
    a₀ = (a₁ * a₂ * a₃)^4 * inv(b₀, A)
    X = a₄ * a₅ * a₃ * a₄ # from Primer
    b₁ = inv(X, A) * a₀ * X # from Primer
    Y = a₂ * a₃ * a₁ * a₂
    return inv(Y, A) * b₁ * Y # b₂ from Primer
end

function Ta(λ::Groups.ΡΛ, i::Integer)
    @assert λ.id == :λ
    return λ[mod1(λ.N - 2i + 1, λ.N), mod1(λ.N - 2i + 2, λ.N)]
end

function Tα(λ::Groups.ΡΛ, i::Integer)
    @assert λ.id == :λ
    return inv(λ[mod1(λ.N - 2i + 2, λ.N), mod1(λ.N - 2i + 1, λ.N)], λ.A)
end

function Te(λ::ΡΛ, ϱ::ΡΛ, i, j)
    @assert i ≠ j

    @assert λ.N == ϱ.N
    @assert λ.A == ϱ.A
    @assert λ.id == :λ && ϱ.id == :ϱ

    @assert iseven(λ.N)
    genus = λ.N÷2
    i = mod1(i, genus)
    j = mod1(j, genus)

    @assert 1 ≤ i ≤ λ.N
    @assert 1 ≤ j ≤ λ.N

    A = λ.A

    if mod(j - (i + 1), genus) == 0
        return Te_diagonal(λ, ϱ, i)
    else
        return inv(Te_lantern(
                A,
                # Our notation:               # Primer notation:
                inv(Ta(λ, i + 1), A),         # b₀
                inv(Ta(λ, i), A),             # a₁
                inv(Tα(λ, i), A),             # a₂
                inv(Te_diagonal(λ, ϱ, i), A), # a₃
                inv(Tα(λ, i + 1), A),         # a₄
                inv(Te(λ, ϱ, i + 1, j), A),   # a₅
            ), A)
    end
end

"""
    rotation_element(G::AutomorphismGroup{<:FreeGroup})
Return the element of `G` which corresponds to shifting generators of the free group.

In the corresponding mapping class group this element acts by rotation of the surface anti-clockwise.
"""
function rotation_element(G::AutomorphismGroup{<:FreeGroup})

    A = alphabet(G)
    @assert iseven(ngens(object(G)))
    genus = ngens(object(G)) ÷ 2

    λ = ΡΛ(:λ, A, 2genus)
    ϱ = ΡΛ(:ϱ, A, 2genus)

    return G(rotation_element(λ, ϱ))
end

function rotation_element(λ::ΡΛ, ϱ::ΡΛ)
    @assert iseven(λ.N)
    genus = λ.N÷2
    A = λ.A

    halftwists = map(1:genus-1) do i
        j = i + 1
        x = Ta(λ, j) * inv(Ta(λ, i), A) * Tα(λ, j) * Te_diagonal(λ, ϱ, i)
        δ = x * Tα(λ, i) * inv(x, A)
        c =
            inv(Ta(λ, j), A) *
            Te(λ, ϱ, i, j) *
            Tα(λ, i)^2 *
            inv(δ, A) *
            inv(Ta(λ, j), A) *
            Ta(λ, i) *
            δ
        z =
            Te_diagonal(λ, ϱ, i) *
            inv(Ta(λ, i), A) *
            Tα(λ, i) *
            Ta(λ, i) *
            inv(Te_diagonal(λ, ϱ, i), A)

        Ta(λ, i) * inv(Ta(λ, j) * Tα(λ, j), A)^6 * (Ta(λ, j) * Tα(λ, j) * z)^4 * c
    end

    τ = (Ta(λ, 1) * Tα(λ, 1))^6 * prod(halftwists)
    return τ
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

struct SymplecticMappingClass{T, F} <: GSymbol
    id::Symbol # :A, :B
    i::UInt
    j::UInt
    minus::Bool
    inv::Bool
    autFn_word::T
    f::F
end

Base.:(==)(a::SymplecticMappingClass, b::SymplecticMappingClass) = a.autFn_word == b.autFn_word

Base.hash(a::SymplecticMappingClass, h::UInt) = hash(a.autFn_word, h)

function SymplecticMappingClass(
    sautFn::AutomorphismGroup{<:FreeGroup},
    id::Symbol,
    i::Integer,
    j::Integer;
    minus = false,
    inverse = false,
)
    @assert i > 0 && j > 0
    id === :A && @assert i ≠ j
    @assert iseven(ngens(object(sautFn)))
    genus = ngens(object(sautFn))÷2

    A = alphabet(sautFn)
    λ = ΡΛ(:λ, A, 2genus)
    ϱ = ΡΛ(:ϱ, A, 2genus)

    w = if id === :A
        Te(λ, ϱ, i, j) *
        inv(Ta(λ, i), A) *
        Tα(λ, i) *
        Ta(λ, i) *
        inv(Te(λ, ϱ, i, j), A) *
        inv(Tα(λ, i), A) *
        inv(Ta(λ, j), A)
    elseif id === :B
        if !minus
            if i ≠ j
                x = Ta(λ, j) * inv(Ta(λ, i), A) * Tα(λ, j) * Te(λ, ϱ, i, j)
                δ = x * Tα(λ, i) * inv(x, A)
                Tα(λ, i) * Tα(λ, j) * inv(δ, A)
            else
                inv(Tα(λ, i), A)
            end
        else
            if i ≠ j
                Ta(λ, i) * Ta(λ, j) * inv(Te(λ, ϱ, i, j), A)
            else
                Ta(λ, i)
            end
        end
    else
        throw("Type not recognized: $id")
    end

    # w is a word defined in the context of A (= alphabet(sautFn))
    # so this "coercion" is correct
    a = sautFn(w)

    f = compiled(a)
    # f = t -> evaluate!(t, a)

    res = SymplecticMappingClass(id, UInt(i), UInt(j), minus, inverse, a, f)

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
    inv_w = inv(m.autFn_word)
    # f(t) = evaluate!(t, inv_w)
    f = compiled(inv_w)
    return SymplecticMappingClass(m.id, m.i, m.j, m.minus, !m.inv, inv_w, f)
end

function evaluate!(
    t::NTuple{N,T},
    smc::SymplecticMappingClass,
    tmp=nothing,
) where {N,T}
    t = smc.f(t)
    for i in 1:N
        normalform!(t[i])
    end
    return t
end
