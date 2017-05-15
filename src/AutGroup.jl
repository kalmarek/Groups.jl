###############################################################################
#
#   AutSymbol/ AutGroup / AutGroupElem
#
###############################################################################

immutable AutSymbol <: GSymbol
   str::String
   pow::Int
   ex::Expr
   func::Function
end

typealias AutGroupElem GWord{AutSymbol}

type AutGroup <: FPGroup
   objectGroup::Group
   gens::Vector{AutSymbol}
end

export AutSymbol, AutGroupElem, AutGroup

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

elem_type(::AutGroup) = AutGroupElem

parent_type(::AutGroupElem) = AutGroup

###############################################################################
#
#   AutSymbol defining functions
#
###############################################################################

function ϱ(i,j, pow=1)
    # @assert i ≠ j
    return v -> [(k==i ? v[i]*v[j]^pow : v[k]) for k in eachindex(v)]
end

function λ(i,j, pow=1)
    # @assert i ≠ j
    return v -> [(k==i ? v[j]^pow*v[i] : v[k]) for k in eachindex(v)]
end

function σ(p::perm, pow=1)
   return v -> [v[(p^pow)[k]] for k in eachindex(v)]
end

ɛ(i, pow=1) = v -> [(k==i ? v[k]^(-1*(2+pow%2)%2) : v[k]) for k in eachindex(v)]

# taken from ValidatedNumerics, under under the MIT "Expat" License:
# https://github.com/JuliaIntervals/ValidatedNumerics.jl/blob/master/LICENSE.md
function subscriptify(n::Int)
    subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
    return join([Char(subscript_0 + i) for i in reverse(digits(n))])
end

function id_autsymbol()
   return AutSymbol("(id)", 0, :(id()), identity)
end

function rmul_autsymbol(i, j; pow::Int=1)
    str = "ϱ"*subscriptify(i)*subscriptify(j)
    return AutSymbol(str, pow, :(ϱ($i, $j, $pow)), ϱ(i, j, pow))
end

function lmul_autsymbol(i, j; pow::Int=1)
    str = "λ"*subscriptify(i)*subscriptify(j)
    return AutSymbol(str, pow, :(λ($i, $j, $pow)), λ(i, j, pow))
end

function flip_autsymbol(i; pow::Int=1)
    str = "ɛ"*subscriptify(i)
    pow = (2+pow%2)%2
    return AutSymbol(str, pow, :(ɛ($i, $pow)), ɛ(i, pow))
end

function perm_autsymbol(p::perm; pow::Int=1)
    if p == parent(p)()
        return id_autsymbol()
    else
        str = "σ"*join([subscriptify(i) for i in p.d])
        p = p^pow
        return AutSymbol(str, 1, :(σ($(p.d), 1)), σ(p, 1))
    end
end

function getperm(s::AutSymbol)
   if s.ex.args[1] == :σ
      p = s.ex.args[2]
      return PermutationGroup(length(p))(p)
   else
      throw(ArgumentError("$s is not a permutation automorphism!"))
   end
end

###############################################################################
#
#   AutGroup / AutGroupElem constructors
#
###############################################################################

function AutGroup(G::FreeGroup; outer=false, special=false)
   n = length(G.gens)
   n == 0 && return AutGroup(G, AutSymbol[])
   S = AutSymbol[]
   indexing = [[i,j] for i in 1:n for j in 1:n if i≠j]
   rmuls = [rmul_autsymbol(i,j) for (i,j) in indexing]
   append!(S, rmuls)
   lmuls = [lmul_autsymbol(i,j) for (i,j) in indexing]
   append!(S, lmuls)
   if !special
      flips = [flip_autsymbol(i) for i in 1:n]
      append!(S, flips)
   end
   if !outer
      perms = collect(elements(PermutationGroup(n)))
      perms = [perm_autsymbol(p) for p in perms[2:end]] # leave the identity
      append!(S, perms)
   end
   return AutGroup(G, S)
end

###############################################################################
#
#   Types call overloads
#
###############################################################################

function (f::AutSymbol){T}(v::Vector{GWord{T}})
    if f.pow == 0
        return v
    end
    return f.func(v)
end

function (F::AutGroupElem)(v::Vector)
    for f in F.symbols
        v = f(v)
    end
    return v
end

function (G::AutGroup)(f::AutSymbol)
   g = AutGroupElem(f)
   g.parent = G
   return g
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

hash(s::AutSymbol, h::UInt) = hash(s.str, hash(s.pow, hash(:AutSymbol, h)))

function change_pow(s::AutSymbol, n::Int)
    if n == 0
        return id_autsymbol()
    end
    symbol = s.ex.args[1]
    if symbol == :ɛ
        return flip_autsymbol(s.ex.args[2], pow=n)
    elseif symbol == :σ
        return perm_autsymbol(s.ex.args[2], pow=n)
    elseif symbol == :ϱ
        s.ex.args[2:end-1]
        return rmul_autsymbol(s.ex.args[2:end-1]..., pow=n)
    elseif symbol == :λ
        return lmul_autsymbol(s.ex.args[2:end-1]..., pow=n)
    elseif symbol == :id
        return s
    else
        warn("Changing power of an unknown type of symbol! $s")
        return AutSymbol(s.str, n, s.ex, s.func)
    end
end

generators(G::AutGroup) = [G(AutGroupElem(elt)) for elt in G.generators]

###############################################################################
#
#   String I/O
#
###############################################################################

###############################################################################
#
#   Comparison
#
###############################################################################

(==)(s::AutSymbol, t::AutSymbol) = s.str == t.str && s.pow == t.pow

###############################################################################
#
#   Binary operators
#
###############################################################################

###############################################################################
#
#   Inversion
#
###############################################################################

inv(f::AutSymbol) = change_pow(f, -f.pow)

###############################################################################
#
#   Misc
#
###############################################################################

function simplify_perms!(W::AutGroupElem)
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
                W.symbols[i] = id_autsymbol()
                W.symbols[i+1] = symmetric_AutSymbol(array(p1*p2))
            end
        end
    end
    deleteat!(W.symbols, find(x -> x.pow == 0, W.symbols))
    return reduced
end

function reduce!(W::AutGroupElem)
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
