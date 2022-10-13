function gersten_alphabet(n::Integer; commutative::Bool = true)
    indexing = [(i, j) for i in 1:n for j in 1:n if i ≠ j]
    S = [ϱ(i, j) for (i, j) in indexing]

    if !commutative
        append!(S, [λ(i, j) for (i, j) in indexing])
    end

    return Alphabet(S)
end

function _commutation_rule(
    ::Type{W},
    A::Alphabet,
    x::S,
    y::S,
) where {S,T,W<:AbstractWord{T}}
    return W(T[A[x], A[y]]) => W(T[A[y], A[x]])
end

function _pentagonal_rule(
    ::Type{W},
    A::Alphabet,
    x::S,
    y::S,
    z::S,
) where {S,T,W<:AbstractWord{T}}
    # x·y·x⁻¹·y⁻¹ => z, i.e. z·y·x => x·y
    return W(T[A[z], A[y], A[x]]) => W(T[A[x], A[y]])
end
function _hexagonal_rule(
    ::Type{W},
    A::Alphabet,
    x::S,
    y::S,
    z::S,
    w::S,
) where {S,T,W<:AbstractWord{T}}
    # x·y⁻¹·z => z·w⁻¹·x
    return W(T[A[x], A[inv(y)], A[z]]) => W(T[A[z], A[w^-1], A[x]])
end

gersten_relations(n::Integer; commutative) =
    gersten_relations(Word{UInt8}, n, commutative = commutative)

function gersten_relations(::Type{W}, n::Integer; commutative) where {W<:AbstractWord}
    @assert n > 1 "Gersten relations are defined only for n>1, got n=$n"
    A = gersten_alphabet(n, commutative = commutative)
    @assert length(A) <= KnuthBendix._max_alphabet_length(W) "Type $W can not represent words over alphabet with $(length(A)) letters."

    rels = Pair{W,W}[]

    for (i, j, k, l) in Iterators.product(1:n, 1:n, 1:n, 1:n)
        if i ≠ j && k ≠ l && k ≠ i && k ≠ j && l ≠ i
            push!(rels, _commutation_rule(W, A, ϱ(i, j), ϱ(k, l)))

            commutative && continue

            push!(rels, _commutation_rule(W, A, λ(i, j), λ(k, l)))
        end
    end

    if !commutative
        for (i, j, k, l) in Iterators.product(1:n, 1:n, 1:n, 1:n)
            if (i ≠ j && k ≠ l && k ≠ j && l ≠ i)
                push!(rels, _commutation_rule(W, A, ϱ(i, j), λ(k, l)))
                push!(rels, _commutation_rule(W, A, λ(i, j), ϱ(k, l)))
            end
        end
    end

    # pentagonal rule:
    # x*y*inv(x)*inv(y)=>z

    for (i, j, k) in Iterators.product(1:n, 1:n, 1:n)
        if (i ≠ j && k ≠ i && k ≠ j)
            push!(rels, _pentagonal_rule(W, A, ϱ(i, j)^-1, ϱ(j, k)^-1, ϱ(i, k)^-1))
            push!(rels, _pentagonal_rule(W, A, ϱ(i, j)^-1, ϱ(j, k), ϱ(i, k)))

            commutative && continue

            push!(rels, _pentagonal_rule(W, A, ϱ(i, j), λ(j, k), ϱ(i, k)^-1))
            push!(rels, _pentagonal_rule(W, A, ϱ(i, j), λ(j, k)^-1, ϱ(i, k)))

            # the same as above, but with ϱ ↔ λ:
            push!(rels, _pentagonal_rule(W, A, λ(i, j)^-1, λ(j, k)^-1, λ(i, k)^-1))
            push!(rels, _pentagonal_rule(W, A, λ(i, j)^-1, λ(j, k), λ(i, k)))

            push!(rels, _pentagonal_rule(W, A, λ(i, j), ϱ(j, k), λ(i, k)^-1))
            push!(rels, _pentagonal_rule(W, A, λ(i, j), ϱ(j, k)^-1, λ(i, k)))
        end
    end

    if !commutative
        for (i, j) in Iterators.product(1:n, 1:n)
            if i ≠ j
                push!(rels, _hexagonal_rule(W, A, ϱ(i, j), ϱ(j, i), λ(i, j), λ(j, i)))
                w = W([A[ϱ(i, j)], A[ϱ(j, i)^-1], A[λ(i, j)]])
                push!(rels, w^2 => inv(w, A)^2)
            end
        end
    end

    return A, rels
end
