struct SurfaceGroup{T, S, R} <: AbstractFPGroup
    genus::Int
    boundaries::Int
    gens::Vector{T}
    relations::Vector{<:Pair{S,S}}
    rws::R
end

include("symplectic_twists.jl")

genus(S::SurfaceGroup) = S.genus

function Base.show(io::IO, S::SurfaceGroup)
    print(io, "π₁ of the orientable surface of genus $(genus(S))")
    if S.boundaries > 0
        print(io, " with $(S.boundaries) boundary components")
    end
end

function SurfaceGroup(genus::Integer, boundaries::Integer)
    @assert genus > 1

    # The (confluent) rewriting systems comes from
    # S. Hermiller, Rewriting systems for Coxeter groups
    # Journal of Pure and Applied Algebra
    # Volume 92, Issue 2, 7 March 1994, Pages 137-148
    # https://doi.org/10.1016/0022-4049(94)90019-1
    # Note: the notation is "inverted":
    # a_g of the article becomes A_g here.

    ltrs = String[]
    for i in 1:genus
        subscript = join('₀'+d for d in reverse(digits(i)))
        append!(ltrs, ["A" * subscript, "a" * subscript, "B" * subscript, "b" * subscript])
    end
    Al = Alphabet(reverse!(ltrs))

    for i in 1:genus
        subscript = join('₀'+d for d in reverse(digits(i)))
        KnuthBendix.set_inversion!(Al, "a" * subscript, "A" * subscript)
        KnuthBendix.set_inversion!(Al, "b" * subscript, "B" * subscript)
    end

    if boundaries == 0
        word = Int[]

        for i in reverse(1:genus)
            x = 4 * i
            append!(word, [x, x-2, x-1, x-3])
        end
        comms = Word(word)
        word_rels = [ comms => one(comms) ]

        rws = KnuthBendix.RewritingSystem(word_rels, KnuthBendix.RecursivePathOrder(Al))
        KnuthBendix.knuthbendix!(rws)
    elseif boundaries == 1
        S = typeof(one(Word(Int[])))
        word_rels = Pair{S, S}[]
        rws = RewritingSystem(word_rels, KnuthBendix.LenLex(Al))
    else
        throw("Not Implemented")
    end

    F = FreeGroup(alphabet(rws))
    rels = [F(lhs)=>F(rhs) for (lhs,rhs) in word_rels]

    return SurfaceGroup(genus, boundaries, [Al[i] for i in 2:2:length(Al)], rels, rws)
end

rewriting(S::SurfaceGroup) = S.rws
KnuthBendix.alphabet(S::SurfaceGroup) = alphabet(rewriting(S))
relations(S::SurfaceGroup) = S.relations

function symplectic_twists(π₁Σ::SurfaceGroup)
    g = genus(π₁Σ)

    saut = SpecialAutomorphismGroup(FreeGroup(2g), max_rules=1000)

    Aij  = [SymplecticMappingClass(saut, :A, i, j) for i in 1:g for j in 1:g if i≠j]

    Bij  = [SymplecticMappingClass(saut, :B, i, j) for i in 1:g for j in 1:g if i≠j]

    mBij = [SymplecticMappingClass(saut, :B, i, j, minus=true) for i in 1:g for j in 1:g if i≠j]

    Bii  = [SymplecticMappingClass(saut, :B, i, i) for i in 1:g]

    mBii = [SymplecticMappingClass(saut, :B, i, i, minus=true) for i in 1:g]

    return [Aij; Bij; mBij; Bii; mBii]
end

KnuthBendix.alphabet(G::AutomorphismGroup{<:SurfaceGroup}) = rewriting(G)

function AutomorphismGroup(π₁Σ::SurfaceGroup; kwargs...)
    S = vcat(symplectic_twists(π₁Σ)...)
    A = Alphabet(S)

    # this is to fix the definitions of symplectic twists:
    # with i->gens(π₁Σ, i) the corresponding automorphisms return
    # reversed words
    domain = ntuple(i->inv(gens(π₁Σ, i)), 2genus(π₁Σ))
    return AutomorphismGroup(π₁Σ, S, A, domain)
end
