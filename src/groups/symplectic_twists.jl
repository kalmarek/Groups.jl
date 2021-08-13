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
        return inv(A, τ) * Te_diagonal(λ, ϱ, 1) * τ
    end

    @assert 1 <= i < n

    NI = (2n - 2i) + 1
    NJ = (2n - 2j) + 1
    I = (2n - 2i) + 2
    J = (2n - 2j) + 2

    g = one(Word(Int[]))
    g *= λ[NJ, NI]                    # β ↦ α*β
    g *= λ[NI, I] * inv(A, ϱ[NI, J])  # α ↦ a*α*b^-1
    g *= inv(A, λ[NJ, NI])            # β ↦ b*α^-1*a^-1*α*β
    g *= λ[J, NI] * inv(A, λ[J, I])   # b ↦ α
    g *= inv(A, λ[J, NI])             # b ↦ b*α^-1*a^-1*α
    g *= inv(A, ϱ[J, NI]) * ϱ[J, I]   # b ↦ b*α^-1*a^-1*α*b*α^-1
    g *= ϱ[J, NI]                     # b ↦ b*α^-1*a^-1*α*b*α^-1*a*α*b^-1

    return g
end

function Te_lantern(A::Alphabet, b₀::T, a₁::T, a₂::T, a₃::T, a₄::T, a₅::T) where {T}
    a₀ = (a₁ * a₂ * a₃)^4 * inv(A, b₀)
    X = a₄ * a₅ * a₃ * a₄ # from Primer
    b₁ = inv(A, X) * a₀ * X # from Primer
    Y = a₂ * a₃ * a₁ * a₂
    return inv(A, Y) * b₁ * Y # b₂ from Primer
end

function Ta(λ::Groups.ΡΛ, i::Integer)
    @assert λ.id == :λ;
    return λ[mod1(λ.N-2i+1, λ.N), mod1(λ.N-2i+2, λ.N)]
end

function Tα(λ::Groups.ΡΛ, i::Integer)
    @assert λ.id == :λ;
    return inv(λ.A, λ[mod1(λ.N-2i+2, λ.N), mod1(λ.N-2i+1, λ.N)])
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
        return inv(A, Te_lantern(
            A,
            # Our notation:               # Primer notation:
            inv(A, Ta(λ, i + 1)),         # b₀
            inv(A, Ta(λ, i)),             # a₁
            inv(A, Tα(λ, i)),             # a₂
            inv(A, Te_diagonal(λ, ϱ, i)), # a₃
            inv(A, Tα(λ, i + 1)),         # a₄
            inv(A, Te(λ, ϱ, i + 1, j)),   # a₅
        ))
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
        x = Ta(λ, j) * inv(A, Ta(λ, i)) * Tα(λ, j) * Te_diagonal(λ, ϱ, i)
        δ = x * Tα(λ, i) * inv(A, x)
        c =
            inv(A, Ta(λ, j)) *
            Te(λ, ϱ, i, j) *
            Tα(λ, i)^2 *
            inv(A, δ) *
            inv(A, Ta(λ, j)) *
            Ta(λ, i) *
            δ
        z =
            Te_diagonal(λ, ϱ, i) *
            inv(A, Ta(λ, i)) *
            Tα(λ, i) *
            Ta(λ, i) *
            inv(A, Te_diagonal(λ, ϱ, i))

        Ta(λ, i) * inv(A, Ta(λ, j) * Tα(λ, j))^6 * (Ta(λ, j) * Tα(λ, j) * z)^4 * c
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
    perm::Vector{Int}
    invperm::Vector{Int}
    f::F
end

Base.:(==)(a::SymplecticMappingClass, b::SymplecticMappingClass) = a.autFn_word == b.autFn_word

Base.hash(a::SymplecticMappingClass, h::UInt) = hash(a.autFn_word, h)

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

    a = sautFn(w)
    g = genus(Σ)
    perm = [2g:-2:1; (2g-1):-2:1]

    f = compiled(a)
    # f = t -> evaluate!(t, a)

    res = SymplecticMappingClass(id, UInt(i), UInt(j), minus, inverse, a, perm, invperm(perm), f)

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
    return SymplecticMappingClass(m.id, m.i, m.j, m.minus, !m.inv, inv_w, m.perm, m.invperm, f)
end

function evaluate!(
    t::NTuple{N,T},
    smc::SymplecticMappingClass,
    tmp=nothing,
) where {N,T}

    t = smc.f(t[smc.perm])[smc.invperm]
    return t
end

struct PermLeftAut <: GSymbol
    perm::Vector{UInt8}

    function PermLeftAut(p::AbstractVector{<:Integer})
        @assert sort(p) == 1:length(p)
        return new(p)
    end
end

function Base.show(io::IO, p::PermLeftAut)
    print(io, "σˡ")
    join(io, (subscriptify(Int(i)) for i in p.perm))
end

Base.inv(p::PermLeftAut) = PermLeftAut(invperm(p.perm))

Base.:(==)(p::PermLeftAut, q::PermLeftAut) = p.perm == q.perm
Base.hash(p::PermLeftAut, h::UInt) = hash(p.perm, hash(PermLeftAut, h))

function evaluate!(v::NTuple{T, N}, p::PermLeftAut, ::Alphabet, tmp=one(first(v))) where {T, N}

    G = parent(first(v))
    A = alphabet(G)
    img = gens(G)[p.perm]
    gens_idcs = Dict(A[A[first(word(im))]] => img[i] for (i,im) in enumerate(gens(G)))

    for elt in v
        copyto!(tmp, elt)
        resize!(word(elt), 0)
        for idx in word(tmp)
            # @show idx
            if haskey(gens_idcs, idx)
                k = gens_idcs[idx]
                append!(word(elt), word(k))
            elseif haskey(gens_idcs, inv(A, idx))
                k = inv(gens_idcs[inv(A, idx)])
                append!(word(elt), word(k))
            else
                push!(word(elt), idx)
            end
        end
        _setnormalform!(elt, false)
        _setvalidhash!(elt, false)

        normalform!(tmp, elt)
        copyto!(elt, tmp)
    end

    return v
end
