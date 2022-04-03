include("eltary_symplectic.jl")

struct SymplecticGroup{N, T, R, A, S} <: MatrixGroup{N,T}
    base_ring::R
    alphabet::A
    gens::S

    function SymplecticGroup{N}(base_ring) where N
        S = symplectic_gens(N, eltype(base_ring))
        alphabet = Alphabet(S)

        return new{
            N,
            eltype(base_ring),
            typeof(base_ring),
            typeof(alphabet),
            typeof(S)
        }(base_ring, alphabet, S)
    end
end

GroupsCore.ngens(Sp::SymplecticGroup) = length(Sp.gens)

Base.show(io::IO, ::SymplecticGroup{N}) where N = print(io, "group of $N×$N symplectic matrices")

function Base.show(
    io::IO,
    ::MIME"text/plain",
    sp::Groups.AbstractFPGroupElement{<:SymplecticGroup{N}}
) where {N}
    Groups.normalform!(sp)
    print(io, "$N×$N symplectic matrix: ")
    KnuthBendix.print_repr(io, word(sp), alphabet(sp))
    println(io)
    Base.print_array(io, matrix_repr(sp))
end

_offdiag_idcs(n) = ((i,j) for i in 1:n for j in 1:n if i ≠ j)

function symplectic_gens(N, T=Int8)
    iseven(N) || throw(ArgumentError("N needs to be even!"))
    n = N÷2

    a_ijs = [ElementarySymplectic{N}(:A, i,j, one(T)) for (i,j) in _offdiag_idcs(n)]
    b_is =  [ElementarySymplectic{N}(:B, n+i,i, one(T)) for i in 1:n]
    c_ijs = [ElementarySymplectic{N}(:B, n+i,j, one(T)) for (i,j) in _offdiag_idcs(n)]

    S = [a_ijs; b_is; c_ijs]

    S = [S; transpose.(S)]

    return unique(S)
end

function _std_symplectic_form(m::AbstractMatrix)
    r,c = size(m)
    r == c || return false
    iseven(r) || return false

    n = r÷2
    𝕆 = zeros(eltype(m), n, n)
    𝕀 = one(eltype(m))*LinearAlgebra.I
    Ω = [𝕆 -𝕀
         𝕀  𝕆]
    return Ω
end

function issymplectic(mat::M, Ω = _std_symplectic_form(mat)) where M <: AbstractMatrix
    return Ω == transpose(mat) * Ω * mat
end
