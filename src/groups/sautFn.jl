include("transvections.jl")
include("gersten_relations.jl")

function SpecialAutomorphismGroup(F::FreeGroup; ordering = KnuthBendix.LenLex, kwargs...)

    n = length(KnuthBendix.alphabet(F)) ÷ 2
    A, rels = gersten_relations(n, commutative = false)
    S = KnuthBendix.letters(A)[1:2(n^2-n)]

    rws = KnuthBendix.RewritingSystem(rels, ordering(A))
    KnuthBendix.knuthbendix!(rws; kwargs...)
    return AutomorphismGroup(F, S, rws, ntuple(i -> gens(F, i), n))
end

KnuthBendix.alphabet(G::AutomorphismGroup{<:FreeGroup}) = alphabet(rewriting(G))

function relations(G::AutomorphismGroup{<:FreeGroup})
    n = length(KnuthBendix.alphabet(object(G))) ÷ 2
    return last(gersten_relations(n, commutative = false))
end

evaluate(f::FPGroupElement{<:AutomorphismGroup{<:FreeGroup}}) = evaluate!(domain(f), f)

function evaluate!(
    t::NTuple{N,T},
    f::FPGroupElement{<:AutomorphismGroup{<:FreeGroup}},
    tmp = one(first(t)),
) where {T<:FPGroupElement,N}
    A = alphabet(f)
    for idx in word(f)
        t = @inbounds evaluate!(t, A[idx], alphabet(object(parent(f))), tmp)::NTuple{N,T}
    end
    return t
end
