include("eltary_matrices.jl")

struct SpecialLinearGroup{N,T,R,A,S} <: MatrixGroup{N,T}
    base_ring::R
    alphabet::A
    gens::S

    function SpecialLinearGroup{N}(base_ring) where {N}
        S = [ElementaryMatrix{N}(i, j, one(base_ring)) for i in 1:N for j in 1:N if i ≠ j]
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

GroupsCore.ngens(SL::SpecialLinearGroup{N}) where {N} = N^2 - N

Base.show(io::IO, SL::SpecialLinearGroup{N,T}) where {N,T} =
    print(io, "special linear group of $N×$N matrices over $T")

function Base.show(
    io::IO,
    ::MIME"text/plain",
    sl::Groups.AbstractFPGroupElement{<:SpecialLinearGroup{N}}
) where {N}

    Groups.normalform!(sl)

    print(io, "SL{$N,$(eltype(sl))} matrix: ")
    KnuthBendix.print_repr(io, word(sl), alphabet(sl))
    println(io)
    Base.print_array(io, matrix_repr(sl))
end

Base.show(io::IO, sl::Groups.AbstractFPGroupElement{<:SpecialLinearGroup}) =
    KnuthBendix.print_repr(io, word(sl), alphabet(sl))
