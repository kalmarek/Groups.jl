###############################################################################
#
#   Naive reduction
#

function freereduce!(::Type{Bool}, w::GWord)
    reduced = true
    for i in 1:syllablelength(w)-1
        s, ns = syllables(w)[i], syllables(w)[i+1]
        if isone(s)
            continue
        elseif s.id == ns.id
            reduced = false
            setmodified!(w)
            p1 = s.pow
            p2 = ns.pow

            syllables(w)[i+1] = change_pow(s, p1 + p2)
            syllables(w)[i] = change_pow(s, 0)
        end
    end
    filter!(!isone, syllables(w))
    return reduced
end

function freereduce!(w::GWord)
    reduced = false
    while !reduced
        reduced = freereduce!(Bool, w)
    end
    return w
end

reduce!(w::GWord) = freereduce!(w)

@doc doc"""
    reduce(w::GWord)
> performs reduction/simplification of a group element (word in generators).
> The default reduction is the free group reduction
> More specific procedures should be dispatched on `GWord`s type parameter.

"""
reduce(w::GWord) = reduce!(deepcopy(w))
