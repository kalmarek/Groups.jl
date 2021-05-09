function KnuthBendix.Alphabet(S::AbstractVector{<:GSymbol})
    S = unique!([S; inv.(S)])
    inversions = [findfirst(==(inv(s)), S) for s in S]
    return Alphabet(S, inversions)
end

struct AutomorphismGroup{G<:Group, T, R, S} <: AbstractFPGroup
    group::G
    gens::Vector{T}
    rws::R
    domain::S
end

object(G::AutomorphismGroup) = G.group

function SpecialAutomorphismGroup(F::FreeGroup;
    ordering=KnuthBendix.LenLex, kwargs...)

    n = length(KnuthBendix.alphabet(F))÷2
    A, rels = gersten_relations(n, commutative=false)
    S = KnuthBendix.letters(A)[1:2(n^2 - n)]

    rws = KnuthBendix.RewritingSystem(rels, ordering(A))
    KnuthBendix.knuthbendix!(rws; kwargs...)
    return AutomorphismGroup(F, S, rws, ntuple(i->gens(F, i), n))
end

KnuthBendix.alphabet(G::AutomorphismGroup{<:FreeGroup}) = alphabet(rewriting(G))
rewriting(G::AutomorphismGroup) = G.rws

function relations(G::AutomorphismGroup)
    n = length(KnuthBendix.alphabet(object(G)))÷2
    return last(gersten_relations(n, commutative=false))
end

_hashing_data(f::FPGroupElement{<:AutomorphismGroup}) = normalform!.(evaluate(f))

function Base.:(==)(g::A, h::A) where A<:FPGroupElement{<:AutomorphismGroup}
    @assert parent(g) === parent(h)

    if _isvalidhash(g) && _isvalidhash(h)
        hash(g) != hash(h) && return false
    end

    normalform!(g)
    normalform!(h)
    word(g) == word(h) && return true

    @assert isnormalform(g)
    @assert isnormalform(h)

    img_computed, imh_computed = false, false

    if !_isvalidhash(g)
        img = _hashing_data(g)
        _update_savedhash!(g, img)
        img_computed = true
    end
    if !_isvalidhash(h)
        imh = _hashing_data(h)
        _update_savedhash!(h, imh)
        imh_computed = true
    end

    @assert _isvalidhash(g)
    @assert _isvalidhash(h)

    hash(g) != hash(h) && return false

    # words are different, but hashes agree
    if !img_computed
        img = _hashing_data(g)
    end
    if !imh_computed
        imh = _hashing_data(h)
    end

    res = img == imh
    !res && @warn "hash collision in == :" g h

    return res
end

# eye-candy

Base.show(io::IO, ::Type{<:FPGroupElement{<:AutomorphismGroup{T}}}) where T <: FreeGroup = print(io, "Automorphism{$T}")

## Automorphism Evaluation

domain(f::FPGroupElement{<:AutomorphismGroup}) = deepcopy(parent(f).domain)
# tuple(gens(object(parent(f)))...)

evaluate(f::FPGroupElement{<:AutomorphismGroup{<:FreeGroup}}) =
    evaluate!(domain(f), f)

function evaluate!(t::NTuple{N, T}, f::FPGroupElement{<:AutomorphismGroup{<:FreeGroup}}) where {T<:FPGroupElement, N}
    A = alphabet(f)
    for idx in word(f)
        t = evaluate!(t, A[idx])::NTuple{N, T}
    end
    return t
end

function evaluate!(v::NTuple{N, T}, s::AutSymbol) where {N, T}
    @assert s.pow in (-1, 1)
    return evaluate!(v, s.fn, isone(s.pow))::NTuple{N, T}
end

function evaluate!(v, ϱ::RTransvect, flag)
    if flag
        append!(New.word(v[ϱ.i]), New.word(v[ϱ.j]   ))
    else
        append!(New.word(v[ϱ.i]), New.word(v[ϱ.j]^-1))
    end
    _setnormalform!(v[ϱ.i], false)
    _setvalidhash!(v[ϱ.i], false)
    return v
end

function evaluate!(v, λ::LTransvect, flag)
    if flag
        prepend!(New.word(v[λ.i]), New.word(v[λ.j]   ))
    else
        prepend!(New.word(v[λ.i]), New.word(v[λ.j]^-1))
    end
    _setnormalform!(v[λ.i], false)
    _setvalidhash!(v[λ.i], false)
    return v
end