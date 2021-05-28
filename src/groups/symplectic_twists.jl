struct SymplecticMappingClass{N, T} <: GSymbol
    id::Symbol # :A, :B
    i::UInt
    j::UInt
    minus::Bool
    inv::Bool
    images::NTuple{N, T}
    invimages::NTuple{N, T}

    function SymplecticMappingClass{N}(G, id, i, j, minus=false, inv=false) where N
        @assert i > 0 && j > 0
        id === :A && @assert i ≠ j

        g = if id === :A
            Te(G, i, j) *
            Ta(N, i)^-1 *
            Tα(N, i) *
            Ta(N, i) *
            Te(G, i, j)^-1 *
            Tα(N,i)^-1 *
            Ta(N, j)^-1
        elseif id === :B
            if !minus
                if i ≠ j
                    x = Ta(N, j) * Ta(N, i)^-1 * Tα(N, j) * Te(G,i,j)
                    δ = x * Tα(N, i) * x^-1
                    Tα(N, i) * Tα(N, j) * inv(δ)
                else
                    Tα(N, i)^-1
                end
            else
                if i ≠ j
                    Ta(N, i) * Ta(N, j) * Te(G, i, j)^-1
                else
                    Ta(N, i)
                end
            end
        else
            throw("Type not recognized: $id")
        end

        res = new(id, i, j, minus, inv,


        )

        return res
    end
end

_indexing(n) = [(i, j) for i = 1:n for j in 1:n if i ≠ j]
_indexing_increasing(n) = [(i, j) for i = 1:n for j = i+1:n]

_λs(N, A) = [ (i == j ? "aaaarggh..." : Word([A[λ(i, j)]])) for i = 1:N, j = 1:N]
_ϱs(N, A) = [ (i == j ? "aaaarggh..." : Word([A[ϱ(i, j)]])) for i = 1:N, j = 1:N]

function Te_diagonal(G, i::Integer)
    N = ngens(object(G))
    # @assert N == size(λ, 1) == size(ϱ, 1)
    @assert iseven(N)
    n = N ÷ 2
    j = i + 1
    @assert 1 <= i < n

    A = KnuthBendix.alphabet(G)
    λ = _λs(N, A)
    ϱ = _ϱs(N, A)

    # comments are for i,j = 1,2
    g = one(word_type(G))
    g *= λ[n+j, n+i]                   # β ↦ α*β
    g *= λ[n+i, i] * inv(A, ϱ[n+i, j]) # α ↦ a*α*b^-1
    g *= inv(A, λ[n+j, n+i])           # β ↦ b*α^-1*a^-1*α*β
    g *= λ[j, n+i] * inv(A, λ[j, i])   # b ↦ α
    g *= inv(A, λ[j, n+i])             # b ↦ b*α^-1*a^-1*α
    g *= inv(A, ϱ[j, n+i]) * ϱ[j, i]   # b ↦ b*α^-1*a^-1*α*b*α^-1
    g *= ϱ[j, n+i]                     # b ↦ b*α^-1*a^-1*α*b*α^-1*a*α*b^-1
    return G(g)
end

function Te_lantern(b₀::T, a₁::T, a₂::T, a₃::T, a₄::T, a₅::T) where {T}
    a₀ = (a₁ * a₂ * a₃)^4 * b₀^-1
    X = a₄ * a₅ * a₃ * a₄
    b₁ = X^-1 * a₀ * X
    Y = a₂ * a₃ * a₁ * a₂
    return Y^-1 * b₁ * Y # b₂
end

Ta(N, i::Integer) = λ[N÷2+i, i]
Tα(N, i::Integer, λ, A) = inv(A, λ[i, N÷2+i])

function Te(G, i, j)
    @assert i ≠ j
    i, j = i < j ? (i, j) : (j, i)

    N = ngens(object(G))

    A = KnuthBendix.alphabet(G)
    λ = _λs(N, A)
    ϱ = _ϱs(N, A)

    if j == i + 1
        return Te_diagonal(G, i)
    else
        return Te_lantern(
            Ta(N, i + 1, λ),
            Ta(N, i, λ),
            Tα(N, i, λ, A),
            Te(N, i, i + 1),
            Tα(N, i + 1, λ, A),
            Te(N, i + 1, j),
        )
    end
end

function mcg_twists(genus::Integer)
    genus < 3 && throw("Not Implemented: genus = $genus < 3")

    G = SpecialAutomorphismGroup(FreeGroup(2genus))
    A = KnuthBendix.alphabet(G)

    λ = _λs(G)
    ϱ = _ϱs(G)

    Tas = [Ta(G, i, λ) for i in 1:genus]
    Tαs = [Tα(G, i, λ, A) for i in 1:genus]

    Tes = [Te(G, i, j, λ, ϱ) for (i,j) in _indexing_increasing(genus)]

    return Tas, Tαs, Tes
end
