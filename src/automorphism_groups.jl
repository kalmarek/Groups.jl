using Permutations

import Base: convert
export AutSymbol, AutWord, rmul_AutSymbol, lmul_AutSymbol, flip_AutSymbol, symmetric_AutSymbol

immutable AutSymbol <: GSymbol
    gen::String
    pow::Int
    ex::Expr
end

(==)(s::AutSymbol, t::AutSymbol) = s.gen == t.gen && s.pow == t.pow
hash(s::AutSymbol, h::UInt) = hash(s.gen, hash(s.pow, hash(:AutSymbol, h)))

IdSymbol(::Type{AutSymbol}) = AutSymbol("(id)", 0, :(IdAutomorphism(N)))

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
    elseif symbol == :IdAutomorphism
        return s
    else
        warn("Changing an unknown type of symbol! $s")
        return AutSymbol(s.gen, n, s.ex)
    end
end

inv(f::AutSymbol) = change_pow(f, -1*f.pow)

function rmul_AutSymbol(i,j; pow::Int=1)
    gen = string('ϱ',Char(8320+i), Char(8320+j)...)
    return AutSymbol(gen, pow, :(ϱ($i,$j)))
end

function lmul_AutSymbol(i,j; pow::Int=1)
    gen = string('λ',Char(8320+i), Char(8320+j)...)
    return AutSymbol(gen, pow, :(λ($i,$j)))
end

function flip_AutSymbol(j; pow::Int=1)
    gen = string('ɛ', Char(8320 + j))
    return AutSymbol(gen, (2+ pow%2)%2, :(ɛ($j)))
end

function symmetric_AutSymbol(perm::Vector{Int}; pow::Int=1)
    # if perm == collect(1:length(perm))
    #     return one(AutSymbol)
    # end
    perm = Permutation(perm)
    ord = order(perm)
    pow = pow % ord
    perm = perm^pow
    gen = string('σ', [Char(8320 + i) for i in array(perm)]...)
    return AutSymbol(gen, 1, :(σ($(array(perm)))))
end

function getperm(s::AutSymbol)
    if s.ex.args[1] == :σ
        return s.ex.args[2]
    else
        throw(ArgumentError("$s is not a permutation automorphism!"))
    end
end

typealias AutWord GWord{AutSymbol}

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
    return reduced
end
