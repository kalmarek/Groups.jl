include("symplectic_twists.jl")

struct SurfaceGroup{T, S, R} <: AbstractFPGroup
    genus::Int
    boundaries::Int
    gens::Vector{T}
    relations::Vector{<:Pair{S,S}}
    rws::R
end

function SurfaceGroup(genus::Integer, boundaries::Integer)
    @assert genus > 1

    ltrs = String[]
    for i in 1:genus
        subscript = join('₀'+d for d in reverse(digits(i)))
        append!(ltrs, ["a" * subscript, "A" * subscript, "b" * subscript, "B" * subscript])
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
        rels = [ comms => one(comms) ]

        rws = RewritingSystem(rels, KnuthBendix.RecursivePathOrder(Al))
        KnuthBendix.knuthbendix!(rws)
    elseif boundaries == 1
        S = typeof(one(Word(Int[])))
        rels = Pair{S, S}[]
        rws = RewritingSystem(rels, KnuthBendix.LenLex(Al))
    else
        throw("Not Implemented")
    end

    return SurfaceGroup(genus, boundaries, KnuthBendix.letters(Al)[2:2:end], rels, rws)
end

rewriting(G::SurfaceGroup) = G.rws
KnuthBendix.alphabet(G::SurfaceGroup) = alphabet(rewriting(G))
relations(G::SurfaceGroup) = G.relations






function mapping_class_group(genus::Integer, punctures::Integer)
    Σ = surface_group(genus, punctures)



    return New.AutomorphismGroup(Σ, S, rws, ntuple(i -> gens(F, i), n))
end

KnuthBendix.alphabet(G::AutomorphismGroup{<:SurfaceGroup}) = alphabet(rewriting(G))
