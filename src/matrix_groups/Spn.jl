include("eltary_symplectic.jl")

struct SymplecticGroup{N,T,R,S} <: AbstractMatrixGroup{N,T}
    base_ring::R
    alphabet::Alphabet{S}
    gens::Vector{S}

    function SymplecticGroup{N}(base_ring) where {N}
        S = symplectic_gens(N, eltype(base_ring))
        alphabet = Alphabet(S)

        T = eltype(base_ring)
        R = typeof(base_ring)
        St = eltype(S)

        return new{N,T,R,St}(base_ring, alphabet, S)
    end
end

GroupsCore.ngens(Sp::SymplecticGroup) = length(Sp.gens)

Base.show(io::IO, ::SymplecticGroup{N,T}) where {N,T} = print(io, "Sp{$N,$T}")

function Base.show(io::IO, ::MIME"text/plain", ::SymplecticGroup{N}) where {N}
    return print(io, "group of $N×$N symplectic matrices")
end

function symplectic_gens(N, T = Int8)
    iseven(N) || throw(ArgumentError("N needs to be even!"))
    n = N ÷ 2

    _offdiag_idcs(n) = ((i, j) for i in 1:n for j in 1:n if i ≠ j)

    a_ijs = [
        ElementarySymplectic{N}(:A, i, j, one(T)) for (i, j) in _offdiag_idcs(n)
    ]
    b_is = [ElementarySymplectic{N}(:B, n + i, i, one(T)) for i in 1:n]
    c_ijs = [
        ElementarySymplectic{N}(:B, n + i, j, one(T)) for
        (i, j) in _offdiag_idcs(n)
    ]

    S = [a_ijs; b_is; c_ijs]

    S = [S; transpose.(S)]

    return unique(S)
end

function _std_symplectic_form(m::AbstractMatrix)
    r, c = size(m)
    r == c || return false
    iseven(r) || return false

    n = r ÷ 2
    𝕆 = zeros(eltype(m), n, n)
    𝕀 = one(eltype(m)) * LinearAlgebra.I
    Ω = [
        𝕆 -𝕀
        𝕀 𝕆
    ]
    return Ω
end

function issymplectic(
    mat::M,
    Ω = _std_symplectic_form(mat),
) where {M<:AbstractMatrix}
    return Ω == transpose(mat) * Ω * mat
end
