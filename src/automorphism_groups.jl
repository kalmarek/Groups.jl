using Permutations

import Base: convert
export AutSymbol, AutWord, rmul_AutSymbol, lmul_AutSymbol, flip_AutSymbol, symmetric_AutSymbol

immutable AutSymbol <: GSymbol
    gen::String
    pow::Int
    ex::Expr
    fmap::Function
    imap::Function
end

function (f::AutSymbol){T}(v::Vector{GWord{T}})
    if f.pow > 0
        map = f.fmap
    else
        map = f.imap
    end
    for i in 1:abs(f.pow)
         v::Vector{GWord{T}} = map(v)
    end
    return v
end

(==)(s::AutSymbol, t::AutSymbol) = s.gen == t.gen && s.pow == t.pow
hash(s::AutSymbol, h::UInt) = hash(s.gen, hash(s.pow, hash(:AutSymbol, h)))

IdSymbol(::Type{AutSymbol}) = AutSymbol("(id)", 0, :(Id(N)), v -> Vector{GWord}(v), v -> Vector{GWord}(v))

function change_pow(s::AutSymbol, n::Int)
    if n == 0
        return one(s)
    end
    symbol = s.ex.args[1]
    if symbol == :ɛ
        return flip_AutSymbol(s.ex.args[2], pow=n)
    elseif symbol == :σ
        return symmetric_AutSymbol(s.ex.args[2], pow=n)
    elseif symbol == :ϱ
        return rmul_AutSymbol(s.ex.args[2], s.ex.args[3], pow=n)
    elseif symbol == :λ
        return lmul_AutSymbol(s.ex.args[2], s.ex.args[3], pow=n)
    elseif symbol == :Id
        return s
    else
        warn("Changing an unknown type of symbol! $s")
        return AutSymbol(s.gen, n, s.ex, s.fmap, s.imap)
    end
end

inv(f::AutSymbol) = change_pow(f, -f.pow)

function ϱ(i,j)
    # @assert i ≠ j
    return v -> [(k!=i ? GWord(v[k]) : v[i]*v[j]) for k in eachindex(v)]
end

function ϱ_inv(i,j)
    # @assert i ≠ j
    return v -> [(k!=i ? GWord(v[k]) : v[i]*v[j]^-1) for k in eachindex(v)]
end


function λ(i,j)
    # @assert i ≠ j
    return v -> ([(k!=i ? GWord(v[k]) : v[j]*v[i]) for k in eachindex(v)])
end

function λ_inv(i,j)
    # @assert i ≠ j
    return v -> ([(k!=i ? GWord(v[k]) : v[j]^-1*v[i]) for k in eachindex(v)])
end


ɛ(i) = v -> [(k!=i ? GWord(v[k]) : v[k]^-1) for k in eachindex(v)]

function σ(perm)
    # @assert sort(perm) == collect(1:length(perm))
    return v -> [GWord(v[perm[k]]) for k in eachindex(v)]
end

function rmul_AutSymbol(i,j; pow::Int=1)
    gen = string('ϱ',Char(8320+i), Char(8320+j)...)
    return AutSymbol(gen, pow, :(ϱ($i,$j)), ϱ(i,j), ϱ_inv(i,j))
end

function lmul_AutSymbol(i,j; pow::Int=1)
    gen = string('λ',Char(8320+i), Char(8320+j)...)
    return AutSymbol(gen, pow, :(λ($i,$j)), λ(i,j), λ_inv(i,j))
end

function flip_AutSymbol(j; pow::Int=1)
    gen = string('ɛ', Char(8320 + j))
    return AutSymbol(gen, (2+ pow%2)%2, :(ɛ($j)), ɛ(j), ɛ(j))
end

function symmetric_AutSymbol(perm::Vector{Int}; pow::Int=1)
    perm = Permutation(perm)
    ord = order(perm)
    pow = pow % ord
    perm = perm^pow
    p = array(perm)
    if p == collect(1:length(p))
        return one(AutSymbol)
    else
        gen = string('σ', [Char(8320 + i) for i in p]...)
        return AutSymbol(gen, 1, :(σ($p)), σ(p), σ(array(inv(perm))))
    end
end

function getperm(s::AutSymbol)
    if s.ex.args[1] == :σ
        return s.ex.args[2]
    else
        throw(ArgumentError("$s is not a permutation automorphism!"))
    end
end

typealias AutWord GWord{AutSymbol}

function (F::AutWord)(v)
    for f in F.symbols
        v = f(v)
    end
    return v
end

convert(::Type{AutWord}, s::AutSymbol) = GWord(s)

function simplify_perms!(W::AutWord)
    reduced = true
    for i in 1:length(W.symbols) - 1
        current = W.symbols[i]
        if current.ex.args[1] == :σ
            if current.pow != 1
                current = symmetric_AutSymbol(perm(current), pow=current.pow)
            end
            next_s = W.symbols[i+1]
            if  next_s.ex.args[1] == :σ
                reduced = false
                if next_s.pow != 1
                    next_s = symmetric_AutSymbol(perm(next_s), pow=next_s.pow)
                end
                p1 = Permutation(getperm(current))
                p2 = Permutation(getperm(next_s))
                W.symbols[i] = one(AutSymbol)
                W.symbols[i+1] = symmetric_AutSymbol(array(p1*p2))
            end
        end
    end
    deleteat!(W.symbols, find(x -> x.pow == 0, W.symbols))
    return reduced
end

function reduce!(W::AutWord)
    if length(W) < 2
        deleteat!(W.symbols, find(x -> x.pow == 0, W.symbols))
    else
        reduced = false
        while !reduced
            reduced = simplify_perms!(W)
            reduced = join_free_symbols!(W)
            deleteat!(W.symbols, find(x -> x.pow == 0, W.symbols))
        end
    end

    W.modified = false
    W.savedhash = hash(W.symbols,hash(typeof(W)))
    return W
end
