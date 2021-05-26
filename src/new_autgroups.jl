function KnuthBendix.Alphabet(S::AbstractVector{<:GSymbol})
    S = unique!([S; inv.(S)])
    inversions = [findfirst(==(inv(s)), S) for s in S]
    return Alphabet(S, inversions)
end

struct AutomorphismGroup{G<:Group,T,R,S} <: AbstractFPGroup
    group::G
    gens::Vector{T}
    rws::R
    domain::S
end

object(G::AutomorphismGroup) = G.group

function SpecialAutomorphismGroup(F::FreeGroup; ordering = KnuthBendix.LenLex, kwargs...)

    n = length(KnuthBendix.alphabet(F)) ÷ 2
    A, rels = gersten_relations(n, commutative = false)
    S = KnuthBendix.letters(A)[1:2(n^2-n)]

    rws = KnuthBendix.RewritingSystem(rels, ordering(A))
    KnuthBendix.knuthbendix!(rws; kwargs...)
    return AutomorphismGroup(F, S, rws, ntuple(i -> gens(F, i), n))
end

KnuthBendix.alphabet(G::AutomorphismGroup{<:FreeGroup}) = alphabet(rewriting(G))
rewriting(G::AutomorphismGroup) = G.rws

function relations(G::AutomorphismGroup)
    n = length(KnuthBendix.alphabet(object(G))) ÷ 2
    return last(gersten_relations(n, commutative = false))
end

function equality_data(f::FPGroupElement{<:AutomorphismGroup})
    imf = evaluate(f)
    # return normalform!.(imf)

    tmp = one(first(imf))
    for g in imf
        normalform!(tmp, g)
        copyto!(g, tmp)
    end
    return imf
end

function Base.:(==)(g::A, h::A) where {A<:FPGroupElement{<:AutomorphismGroup}}
    @assert parent(g) === parent(h)

    if _isvalidhash(g) && _isvalidhash(h)
        hash(g) != hash(h) && return false
    end

    length(word(g)) > 8 && normalform!(g)
    length(word(h)) > 8 && normalform!(h)

    word(g) == word(h) && return true

    img_computed, imh_computed = false, false

    if !_isvalidhash(g)
        img = equality_data(g)
        _update_savedhash!(g, img)
        img_computed = true
    end
    if !_isvalidhash(h)
        imh = equality_data(h)
        _update_savedhash!(h, imh)
        imh_computed = true
    end

    @assert _isvalidhash(g)
    @assert _isvalidhash(h)

    hash(g) != hash(h) && return false

    # words are different, but hashes agree
    if !img_computed
        img = equality_data(g)
    end
    if !imh_computed
        imh = equality_data(h)
    end

    equal = img == imh
    equal || @warn "hash collision in == :" g h

    return equal
end

function Base.isone(g::FPGroupElement{<:AutomorphismGroup})
    if length(word(g)) > 8
        normalform!(g)
    end
    return evaluate(g) == parent(g).domain
end

# eye-candy

Base.show(io::IO, ::Type{<:FPGroupElement{<:AutomorphismGroup{T}}}) where {T<:FreeGroup} =
    print(io, "Automorphism{$T,…}")

## Automorphism Evaluation

domain(f::FPGroupElement{<:AutomorphismGroup}) = deepcopy(parent(f).domain)
# tuple(gens(object(parent(f)))...)

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
