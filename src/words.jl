syllablelength(w::GWord) = length(w.symbols)
syllables(w::GWord) = w.symbols
ismodified(w::GWord) = w.modified
setmodified!(w::GWord) = (w.modified = true; w)
unsetmodified!(w::GWord) = (w.modified = false; w)
savehash!(w::GWord, h::UInt) = (w.savedhash = h; w)
savedhash(w::GWord) = w.savedhash
parent(w::GWord) = w.parent
hasparent(w::GWord) = isdefined(w, :parent)
setparent!(w::GWord, G::AbstractFPGroup) = (w.parent = G; w)

Base.isempty(w::GWord) = isempty(syllables(w))
Base.isone(w::GWord) = (freereduce!(w); isempty(w))
Base.one(w::GWord) = one(parent(w))

function Base.iterate(w::GWord, state=(syllable=1, pow=1))
    state.syllable > syllablelength(w) && return nothing
    next = iterate(syllables(w)[state.syllable], state.pow)
    next === nothing && return iterate(w, (syllable=state.syllable+1, pow=1))
    return first(next), (syllable=state.syllable, pow=last(next))
end

Base.eltype(::Type{<:GWord{T}}) where T = T
Base.length(w::GWord) = isempty(w) ? 0 : sum(length, syllables(w))
Base.size(w::GWord) = (length(w),)
Base.lastindex(w::GWord) = length(w)

Base.@propagate_inbounds function Base.getindex(w::GWord, i::Integer)
    csum = 0
    idx = 0
    @boundscheck 0 < i <= length(w) || throw(BoundsError(w, i))
    while csum < i
        idx += 1
        csum += length(syllables(w)[idx])
    end
    return first(syllables(w)[idx])
end

Base.@propagate_inbounds Base.getindex(w::GWord, itr) = [w[i] for i in itr]

# no setindex! for syllable based words

Base.convert(::Type{GW}, s::GSymbol) where GW <: GWord = GW(s)
