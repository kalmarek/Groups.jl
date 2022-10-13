include("transvections.jl")
include("gersten_relations.jl")

function SpecialAutomorphismGroup(F::FreeGroup; ordering = KnuthBendix.LenLex, kwargs...)

    n = length(alphabet(F)) ÷ 2
    A, rels = gersten_relations(n, commutative=false)
    S = [A[i] for i in 1:2:length(A)]

    maxrules = 1000*n

    rws = KnuthBendix.RewritingSystem(rels, ordering(A))
    Logging.with_logger(Logging.NullLogger()) do
        # the rws is not confluent, let's suppress warning about it
        KnuthBendix.knuthbendix!(rws; maxrules=maxrules, kwargs...)
    end
    return AutomorphismGroup(F, S, rws, ntuple(i -> gens(F, i), n))
end

KnuthBendix.alphabet(G::AutomorphismGroup{<:FreeGroup}) = alphabet(rewriting(G))

function relations(G::AutomorphismGroup{<:FreeGroup})
    n = length(alphabet(object(G))) ÷ 2
    return last(gersten_relations(n, commutative = false))
end
