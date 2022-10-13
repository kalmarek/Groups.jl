include("transvections.jl")
include("gersten_relations.jl")

function SpecialAutomorphismGroup(F::FreeGroup; ordering = KnuthBendix.LenLex, kwargs...)

    n = length(alphabet(F)) ÷ 2
    A, rels = gersten_relations(n, commutative=false)
    S = [A[i] for i in 1:2:length(A)]

    max_rules = 1000 * n

    rws = Logging.with_logger(Logging.NullLogger()) do
        rws = KnuthBendix.RewritingSystem(rels, ordering(A))
        # the rws is not confluent, let's suppress warning about it
        KnuthBendix.knuthbendix(rws, KnuthBendix.Settings(; max_rules=max_rules, kwargs...))
    end

    idxA = KnuthBendix.IndexAutomaton(rws)
    return AutomorphismGroup(F, S, idxA, ntuple(i -> gens(F, i), n))
end

function relations(G::AutomorphismGroup{<:FreeGroup})
    n = length(alphabet(object(G))) ÷ 2
    return last(gersten_relations(n, commutative = false))
end
