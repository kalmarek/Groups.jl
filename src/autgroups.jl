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
rewriting(G::AutomorphismGroup) = G.rws

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

function Base.:(==)(g::A, h::A) where {A<:AbstractFPGroupElement{<:AutomorphismGroup}}
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

function Base.isone(g::AbstractFPGroupElement{<:AutomorphismGroup})
    if length(word(g)) > 8
        normalform!(g)
    end
    return evaluate(g) == parent(g).domain
end

# eye-candy

Base.show(io::IO, ::Type{<:FPGroupElement{<:AutomorphismGroup{T}}}) where {T} =
    print(io, "Automorphism{$T, …}")

Base.show(io::IO, A::AutomorphismGroup) = print(io, "automorphism group of ", object(A))

function Base.show(io::IO, ::MIME"text/plain", a::AbstractFPGroupElement{<:AutomorphismGroup})
    println(io, " ┌ $(a):")
    d = domain(a)
    im = evaluate(a)
    for (x, imx) in zip(d, im[1:end-1])
        println(io, " │   $x ↦ $imx")
    end
    println(io, " └   $(last(d)) ↦ $(last(im))")
end

## Automorphism Evaluation

domain(f::AbstractFPGroupElement{<:AutomorphismGroup}) = deepcopy(parent(f).domain)
# tuple(gens(object(parent(f)))...)

evaluate(f::AbstractFPGroupElement{<:AutomorphismGroup}) = evaluate!(domain(f), f)

function evaluate!(
    t::NTuple{N,T},
    f::AbstractFPGroupElement{<:AutomorphismGroup{<:Group}},
    tmp = one(first(t)),
) where {N, T<:FPGroupElement}
    A = alphabet(f)
    for idx in word(f)
        t = @inbounds evaluate!(t, A[idx], tmp)::NTuple{N,T}
    end
    return t
end

evaluate!(t::NTuple{N, T}, s::GSymbol, tmp=nothing) where {N, T} = throw("you need to implement `evaluate!(::$(typeof(t)), ::$(typeof(s)), ::Alphabet, tmp=one(first(t)))`")

# forward evaluate by substitution

struct LettersMap{T, A}
    indices_map::Dict{Int, T}
    A::A
end

function LettersMap(a::FPGroupElement{<:AutomorphismGroup})
    dom = domain(a)
    @assert all(isone ∘ length ∘ word, dom)
    A = alphabet(first(dom))
    first_letters = first.(word.(dom))
    img = evaluate!(dom, a)

    # (dom[i] → img[i] is a map from domain to images)
    # we need a map from alphabet indices → (gens, gens⁻¹) → images
    # here we do it for elements of the domain
    # (trusting it's a set of generators that define a)
    @assert length(dom) == length(img)

    indices_map = Dict(A[A[fl]] => word(im) for (fl, im) in zip(first_letters, img))
    # inverses of generators are dealt lazily in getindex

    return LettersMap(indices_map, A)
end


function Base.getindex(lm::LettersMap, i::Integer)
    # here i is an index of an alphabet
    @boundscheck 1 ≤ i ≤ length(lm.A)

    if !haskey(lm.indices_map, i)
        img = if haskey(lm.indices_map, inv(lm.A, i))
            inv(lm.A, lm.indices_map[inv(lm.A, i)])
        else
            @warn "LetterMap: neither $i nor its inverse has assigned value"
            one(valtype(lm.indices_map))
        end
        lm.indices_map[i] = img
    end
    return lm.indices_map[i]
end

function (a::FPGroupElement{<:AutomorphismGroup})(g::FPGroupElement)
    @assert object(parent(a)) === parent(g)
    img_w = evaluate(word(g), LettersMap(a))
    return parent(g)(img_w)
end

evaluate(w::AbstractWord, lm::LettersMap) = evaluate!(one(w), w, lm)

function evaluate!(res::AbstractWord, w::AbstractWord, lm::LettersMap)
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
            k = findfirst(==(inv(A, l)), first_ltrs)
            if k !== nothing
                push!(args[idx].args, :(inv(d[$k])))
                continue
            end
            throw("Letter $l doesn't seem to be mapped anywhere!")
        end
    end
    locals = Dict{Expr, Symbol}()
    locals_counter = 0
    for (i,v) in enumerate(args)
        @assert length(v.args) >= 2
        if length(v.args) > 2
            for (j, a) in pairs(v.args)
                if a isa Expr &&  a.head == :call "$a"
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
        $([:(local $v = $k) for (k,v) in locals]...)
    end

    # return args, locals

    return :(d -> begin
        @boundscheck @assert length(d) == $(length(d))
        $q
        @inbounds tuple($(args...))
    end)
end
