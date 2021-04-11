###############################################################################
#
#   Replacement of symbols / sub-words
#

issubsymbol(s::GSymbol, t::GSymbol) =
    s.id == t.id && (0 ≤ s.pow ≤ t.pow || 0 ≥ s.pow ≥ t.pow)

function issubsymbol(s::FreeSymbol, w::GWord, sindex::Integer)
    @boundscheck 1 ≤ sindex ≤ syllablelength(w) || throw(BoundsError(w, sindex))
    return issubsymbol(s, syllables(w)[sindex])
end

function issubword(z::GWord, w::GWord, sindex::Integer)
    isempty(z) && return true
    @boundscheck 1 ≤ sindex ≤ syllablelength(w) || throw(BoundsError(w, sindex))
    n = syllablelength(z)
    n == 1 && return issubsymbol(first(syllables(z)), syllables(w)[sindex])

    lastindex = sindex + n - 1
    lastindex > syllablelength(w) && return false

    issubsymbol(first(z), syllables(w)[sindex]) || return false
    issubsymbol(syllables(z)[end], syllables(w)[lastindex]) || return false
    for (zidx, widx) in zip(2:n-1, sindex+1:lastindex-1)
        syllables(z)[zidx] == syllables(w)[widx] || return false
    end
    return true
end

"""

Find the first syllable index k>=i such that Z < syllables(W)[k:k+syllablelength(Z)-1]
"""
function Base.findnext(subword::GWord, word::GWord, start::Integer)
    @boundscheck 1 ≤ start ≤ syllablelength(word) || throw(BoundsError(word, start))
    isempty(subword) && return start
    stop = syllablelength(word) - syllablelength(subword) +1

    for idx in start:1:stop
        issubword(subword, word, idx) && return idx
    end
    return nothing
end

function Base.findnext(s::FreeSymbol, word::GWord, start::Integer)
    @boundscheck 1 ≤ start ≤ syllablelength(word) || throw(BoundsError(word, start))
    isone(s) && return start
    stop = syllablelength(word)

    for idx in start:1:stop
        issubsymbol(s, word, idx) && return idx
    end
    return nothing
end

function Base.findprev(subword::GWord, word::GWord, start::Integer)
    @boundscheck 1 ≤ start ≤ syllablelength(word) || throw(BoundsError(word, start))
    isempty(subword) && return start
    stop = 1

    for idx in start:-1:1
        issubword(subword, word, idx) && return idx
    end
    return nothing
end

function Base.findprev(s::FreeSymbol, word::GWord, start::Integer)
    @boundscheck 1 ≤ start ≤ syllablelength(word) || throw(BoundsError(word, start))
    isone(s) && return start
    stop = 1

    for idx in start:-1:stop
        issubsymbol(s, word, idx) && return idx
    end
    return nothing
end

Base.findfirst(subword::GWord, word::GWord) = findnext(subword, word, 1)
Base.findlast(subword::GWord, word::GWord) =
    findprev(subword, word, syllablelength(word)-syllablelength(subword)+1)

function Base.replace!(out::GW, W::GW, lhs_rhs::Pair{GS, T}; count::Integer=typemax(Int)) where
        {GS<:GSymbol, T<:GWord, GW<:GWord}
    (count == 0 || isempty(W)) && return W
    count < 0 && throw(DomainError(count, "`count` must be non-negative."))

    lhs, rhs = lhs_rhs

    sW = syllables(W)
    sW_idx = 1
    r = something(findnext(lhs, W, sW_idx), 0)

    sout = syllables(out)
    resize!(sout, 0)
    sizehint!(sout, syllablelength(W))

    c = 0

    while !iszero(r)
        append!(sout, view(sW, sW_idx:r-1))
        a, b = divrem(sW[r].pow, lhs.pow)

        if b != 0
            push!(sout, change_pow(sW[r], b))
        end

        append!(sout, repeat(syllables(rhs), a))

        sW_idx = r+1
        sW_idx > syllablelength(W) && break

        r = something(findnext(lhs, W, sW_idx), 0)
        c += 1
        c == count && break
    end
    append!(sout, sW[sW_idx:end])
    return freereduce!(out)
end

function Base.replace!(out::GW, W::GW, lhs_rhs::Pair{T, T}; count::Integer=typemax(Int)) where
    {GW<:GWord, T <: GWord}
    (count == 0 || isempty(W)) && return W
    count < 0 && throw(DomainError(count, "`count` must be non-negative."))

    lhs, rhs = lhs_rhs
    lhs_slen = syllablelength(lhs)
    lhs_slen == 1 && return replace!(out, W, first(syllables(lhs))=>rhs; count=count)

    sW = syllables(W)
    sW_idx = 1
    r = something(findnext(lhs, W, sW_idx), 0)

    sout = syllables(out)
    resize!(sout, 0)
    sizehint!(sout, syllablelength(W))

    c = 0

    while !iszero(r)
        append!(sout, view(sW, sW_idx:r-1))

        exp = sW[r].pow - first(syllables(lhs)).pow
        if exp != 0
            push!(sout, change_pow(sW[r], exp))
        end

        append!(sout, syllables(rhs))

        exp = sW[r+lhs_slen-1].pow - last(syllables(lhs)).pow
        if exp != 0
            push!(sout, change_pow(sW[r+lhs_slen-1], exp))
        end

        sW_idx = r+lhs_slen
        sW_idx > syllablelength(W) && break

        r = something(findnext(lhs, W, sW_idx), 0)
        c += 1
        c == count && break
    end

    # copy the rest
    append!(sout, sW[sW_idx:end])
    return freereduce!(out)
end

function Base.replace(W::GW, lhs_rhs::Pair{T, T}; count::Integer=typemax(Int)) where
    {GW<:GWord, T <: GWord}
    return replace!(one(W), W, lhs_rhs; count=count)
end

function Base.replace(W::GW, subst_dict::Dict{T,T}) where {GW<:GWord, T<:GWord}
    out = W
    for toreplace in reverse!(sort!(collect(keys(subst_dict)), by=length))
        replacement = subst_dict[toreplace]
        if length(toreplace) > length(out)
            continue
        end
        out = replace(out, toreplace=>replacement)
    end
    return out
end
