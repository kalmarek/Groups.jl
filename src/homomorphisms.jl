"""
    Homomorphism(f, G::AbstractFPGroup, H::AbstractFPGroup[, check=true])
Struct representing homomorphism map from `G` to `H` given by map `f`.

To define `h = Homomorphism(f, G, H)` function (or just callable) `f` must
implement method `f(i::Integer, source, target)::AbstractWord` with the
following meaning. Suppose that word `w = Word([i])` consists of a single
letter in the `alphabet` of `source` (usually it means that in `G` it
represents a generator or its inverse). Then `f(i, G, H)` must return the
**word** representing the image in `H` of `G(w)` under the homomorphism.

In more mathematical terms it means that if `h(G(w)) == h`, then
`f(i, G, H) == word(h)`.

Images of both `AbstractWord`s and elements of `G` can be obtained by simply
calling `h(w)`, or `h(g)`.

If `check=true` then the correctness of the definition of `h` will be performed
when creating the homomorphism.

!!! note
    `f(i, G, H)` must be implemented for all letters in the alphabet of `G`,
    not only for those `i` which represent `gens(G)`. Function `f` will be
    evaluated exactly once per letter of `alphabet(G)` and the results will be
    cached.

# Examples
```julia
julia> F₂ = FreeGroup(2)
free group on 2 generators

julia> g,h = gens(F₂)
2-element Vector{FPGroupElement{FreeGroup{Symbol, KnuthBendix.LenLex{Symbol}}, …}}:
 f1
 f2

julia> ℤ² = FPGroup(F₂, [g*h => h*g])
Finitely presented group generated by:
        { f1  f2 },
subject to relations:
 f1*f2 => f2*f1

julia> hom = Groups.Homomorphism(
           (i, G, H) -> Groups.word_type(H)([i]),
           F₂,
           ℤ²
       )
Homomorphism
 from : free group on 2 generators
 to   : ⟨ f1  f2 |
          f1*f2 => f2*f1 ⟩

julia> hom(g*h*inv(g))
f2

julia> hom(g*h*inv(g)) == hom(h)
true
```

"""
struct Homomorphism{Gr1, Gr2, I, W}
    gens_images::Dict{I, W}
    source::Gr1
    target::Gr2

    function Homomorphism(
        f,
        source::AbstractFPGroup,
        target::AbstractFPGroup;
        check=true
    )
        A = alphabet(source)
        dct = Dict(i=>convert(word_type(target), f(i, source, target))
            for i in 1:length(A))
        I = eltype(word_type(source))
        W = word_type(target)
        hom = new{typeof(source), typeof(target), I, W}(dct, source, target)

        if check
            @assert hom(one(source)) == one(target)
            for x in gens(source)

                @assert hom(x^-1) == hom(x)^-1

                for y in gens(source)
                    @assert hom(x*y) == hom(x)*hom(y)
                    @assert hom(x*y)^-1 == hom(y^-1)*hom(x^-1)
                end
            end
            for (lhs, rhs) in relations(source)
                relator = lhs * inv(rhs, alphabet(source))
                im_r = hom.target(hom(relator))
                @assert isone(im_r) "Map does not define a homomorphism: h($relator) = $(im_r) ≠ $(one(target))."
            end
        end
        return hom
    end
end

function (h::Homomorphism)(w::AbstractWord)
    result = one(word_type(h.target)) # Word
    for l in w
        append!(result, h.gens_images[l])
    end
    return result
end

function (h::Homomorphism)(g::AbstractFPGroupElement)
    @assert parent(g) === h.source
    w = h(word(g))
    return h.target(w)
end

Base.show(io::IO, h::Homomorphism) = print(io, "Homomorphism\n from : $(h.source)\n to   : $(h.target)")
