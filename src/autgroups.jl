function KnuthBendix.Alphabet(S::AbstractVector{<:GSymbol})
    S = union(S, inv.(S))
    inversions = [findfirst(==(inv(s)), S) for s in S]
    return Alphabet(S, inversions)
end

mutable struct AutomorphismGroup{G<:Group,T,RW,S} <: AbstractFPGroup
    group::G
    gens::Vector{T}
    rw::RW
    domain::S
end

object(G::AutomorphismGroup) = G.group
rewriting(G::AutomorphismGroup) = G.rw

function equality_data(f::AbstractFPGroupElement{<:AutomorphismGroup})
    imf = evaluate(f)
    # return normalform!.(imf)

    tmp = one(first(imf))
    for g in imf
        normalform!(tmp, g)
        copyto!(g, tmp)
    end
    return imf
end

function Base.:(==)(
    g::A,
    h::A,
) where {A<:AbstractFPGroupElement{<:AutomorphismGroup}}
    @assert parent(g) === parent(h)

    if _isvalidhash(g) && _isvalidhash(h)
        hash(g) != hash(h) && return false
    end

    normalform!(g)
    normalform!(h)

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

function Base.isone(g::AbstractFPGroupElement{<:AutomorphismGroup})
    if length(word(g)) > 8
        normalform!(g)
    end
    return evaluate(g) == parent(g).domain
end

# eye-candy

function Base.show(
    io::IO,
    ::Type{<:FPGroupElement{<:AutomorphismGroup{T}}},
) where {T}
    return print(io, "Automorphism{$T, …}")
end

function Base.show(io::IO, A::AutomorphismGroup)
    return print(io, "automorphism group of ", object(A))
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    a::AbstractFPGroupElement{<:AutomorphismGroup},
)
    println(io, " ┌ $(a):")
    d = domain(a)
    im = evaluate(a)
    for (x, imx) in zip(d, im[1:end-1])
        println(io, " │   $x ↦ $imx")
    end
    return println(io, " └   $(last(d)) ↦ $(last(im))")
end

## Automorphism Evaluation

function domain(f::AbstractFPGroupElement{<:AutomorphismGroup})
    return deepcopy(parent(f).domain)
end
# tuple(gens(object(parent(f)))...)

function evaluate(f::AbstractFPGroupElement{<:AutomorphismGroup})
    return evaluate!(domain(f), f)
end

function evaluate!(
    t::NTuple{N,T},
    f::AbstractFPGroupElement{<:AutomorphismGroup{<:Group}},
    tmp = one(first(t)),
) where {N,T<:FPGroupElement}
    A = alphabet(f)
    for idx in word(f)
        t = @inbounds evaluate!(t, A[idx], tmp)::NTuple{N,T}
    end
    return t
end

function evaluate!(t::NTuple{N,T}, s::GSymbol, tmp = nothing) where {N,T}
    throw(
        "you need to implement `evaluate!(::$(typeof(t)), ::$(typeof(s)), ::Alphabet, tmp=one(first(t)))`",
    )
end

# forward evaluate by substitution

struct LettersMap{W<:AbstractWord,A}
    indices_map::Dict{Int,W}
    A::A
end

function LettersMap(a::FPGroupElement{<:AutomorphismGroup})
    dom = domain(a)
    if all(isone ∘ length ∘ word, dom)
        A = alphabet(first(dom))
        first_letters = first.(word.(dom))
        img = evaluate!(dom, a)

        # (dom[i] → img[i] is a map from domain to images)
        # we need a map from alphabet indices → (gens, gens⁻¹) → images
        # here we do it for elements of the domain
        # (trusting it's a set of generators that define a)
        @assert length(dom) == length(img)

        indices_map =
            Dict(Int(fl) => word(im) for (fl, im) in zip(first_letters, img))
        # inverses of generators are dealt lazily in getindex
    else
        throw("LettersMap is not implemented for non-generators in domain")
    end

    return LettersMap(indices_map, A)
end

function Base.getindex(lm::LettersMap{W}, i::Integer) where {W}
    # here i is an index of an alphabet
    @boundscheck 1 ≤ i ≤ length(lm.A)

    if !haskey(lm.indices_map, i)
        I = inv(i, lm.A)
        if haskey(lm.indices_map, I)
            img = inv(lm.indices_map[I], lm.A)
            lm.indices_map[i] = img
        else
            lm.indices_map[i] = W([i])
            lm.indices_map[I] = W([I])
        end
    end
    return lm.indices_map[i]
end

function (a::FPGroupElement{<:AutomorphismGroup})(g::FPGroupElement)
    @assert object(parent(a)) === parent(g)
    img_w = evaluate(word(g), LettersMap(a))
    return parent(g)(img_w)
end

evaluate(w::AbstractWord, lm::LettersMap) = evaluate!(similar(w), w, lm)

function evaluate!(res::AbstractWord, w::AbstractWord, lm::LettersMap)
    resize!(res, 0)
    for i in w
        append!(res, lm[i])
    end
    return res
end

# compile automorphisms

compiled(a) = eval(generated_evaluate(a))

function generated_evaluate(a::FPGroupElement{<:AutomorphismGroup})
    lm = Groups.LettersMap(a)
    d = Groups.domain(a)
    @assert all(length.(word.(d)) .== 1)
    A = alphabet(first(d))
    first_ltrs = first.(word.(d))

    args = [Expr(:call, :*) for _ in first_ltrs]

    for (idx, letter) in enumerate(first_ltrs)
        for l in lm[letter]
            k = findfirst(==(l), first_ltrs)
            if k !== nothing
                push!(args[idx].args, :(d[$k]))
                continue
            end
            k = findfirst(==(inv(l, A)), first_ltrs)
            if k !== nothing
                push!(args[idx].args, :(inv(d[$k])))
                continue
            end
            throw("Letter $l doesn't seem to be mapped anywhere!")
        end
    end
    locals = Dict{Expr,Symbol}()
    locals_counter = 0
    for (i, v) in enumerate(args)
        @assert length(v.args) >= 2
        if length(v.args) > 2
            for (j, a) in pairs(v.args)
                if a isa Expr && a.head == :call
                    "$a"
                    @assert a.args[1] == :inv
                    if !(a in keys(locals))
                        locals[a] = Symbol("var_#$locals_counter")
                        locals_counter += 1
                    end
                    v.args[j] = locals[a]
                end
            end
        else
            args[i] = v.args[2]
        end
    end

    q = quote
        $([:(local $v = $k) for (k, v) in locals]...)
    end

    # return args, locals

    return :(d -> begin
        @boundscheck @assert length(d) == $(length(d))
        $q
        @inbounds tuple($(args...))
    end)
end
