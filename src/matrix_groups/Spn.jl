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
    print(io, "$N×$N Symplectic matrix: ")
    KnuthBendix.print_repr(io, word(sp), alphabet(sp))
    println(io)
    Base.print_array(io, matrix_repr(sp))
end

function symplectic_gens(N, T=Int8)
    iseven(N) || throw(ArgumentError("N needs to be even!"))
    n = N÷2

    a_ijs = [ElementarySymplectic{N}(:A, i,j, one(T)) for (i,j) in offdiagonal_indexing(n)]
    b_is =  [ElementarySymplectic{N}(:B, n+i,i, one(T)) for i in 1:n]
    c_ijs = [ElementarySymplectic{N}(:B, n+i,j, one(T)) for (i,j) in offdiagonal_indexing(n)]

    S = [a_ijs; b_is; c_ijs]

    S = [S; transpose.(S)]

    return unique(S)
end

function _std_symplectic_form(m::AbstractMatrix)
    r,c = size(m)
    r == c || return false
    iseven(r) || return false

    n = r÷2
    Ω = zero(m)
    for i in 1:n
        Ω[2i-1:2i, 2i-1:2i] .= [0 -1; 1 0]
    end
    return Ω
end

function issymplectic(mat::M, Ω = _std_symplectic_form(mat)) where M <: AbstractMatrix
    r, c = size(mat)
    return Ω == transpose(mat) * Ω * mat
end
