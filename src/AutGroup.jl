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

type AutGroup <: AbstractFPGroup
   objectGroup::Group
   gens::Vector{AutSymbol}
end

export AutGroupElem, AutGroup

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
        p = p^pow
        str = "σ"*join([subscriptify(i) for i in p.d])
        return AutSymbol(str, 1, :(σ($(p.d), 1)), σ(p, 1))
    end
end

function perm_autsymbol(a::Vector{Int})
   G = PermutationGroup(length(a))
   return perm_autsymbol(G(a))
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

function AutGroup(G::FreeGroup; special=false)
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
   return AutGroup(G, S)
end

###############################################################################
#
#   Types call overloads
#
###############################################################################

function (G::AutGroup)()
   id = AutGroupElem(id_autsymbol())
   id.parent = G
   return id
end

function (G::AutGroup)(f::AutSymbol)
   g = AutGroupElem(f)
   g.parent = G
   return g
end

function (G::AutGroup)(g::AutGroupElem)
   g.parent = G
   return g
end

###############################################################################
#
#   Functional call overloads for evaluation of AutSymbol and AutGroupElem
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

###############################################################################
#
#   Basic manipulation
#
###############################################################################

hash(s::AutSymbol, h::UInt) = hash(s.str, hash(s.pow, hash(:AutSymbol, h)))

function hash(g::AutGroupElem, h::UInt)
   gs = gens(parent(g).objectGroup)
   return hash(g(gs), hash(typeof(g), hash(parent(g), h)))
end

function change_pow(s::AutSymbol, n::Int)
    if n == 0
        return id_autsymbol()
    end
    symbol = s.ex.args[1]
    if symbol == :ɛ
        return flip_autsymbol(s.ex.args[2], pow=n)
    elseif symbol == :σ
        G = PermutationGroup(length(s.ex.args[2]))
        return perm_autsymbol(G(s.ex.args[2]), pow=n)
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

length(s::AutSymbol) = abs(s.pow)

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::AutGroup)
   print(io, "Automorphism Group of $(G.objectGroup)\n")
   print(io, "Generated by $(join(G.gens, ","))")
end

###############################################################################
#
#   Comparison
#
###############################################################################

(==)(s::AutSymbol, t::AutSymbol) = s.str == t.str && s.pow == t.pow

function (==)(g::AutGroupElem, h::AutGroupElem)
   parent(g) == parent(h) || return false
   G = parent(g).objectGroup
   gs = gens(G)
   return g(gs) == h(gs)
end

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

ispermauto(s::AutSymbol) = s.ex.args[1] == :σ

function simplify_perms!(W::AutGroupElem)
    reduced = true
    for i in 1:length(W.symbols) - 1
        current = W.symbols[i]
        if ispermauto(current)
            if current.pow != 1
                current = perm_autsymbol(perm(current), pow=current.pow)
            end
            next_s = W.symbols[i+1]
            if  ispermauto(next_s)
                reduced = false
                if next_s.pow != 1
                    next_s = perm_autsymbol(perm(next_s), pow=next_s.pow)
                end
                p1 = getperm(current)
                p2 = getperm(next_s)
                W.symbols[i] = id_autsymbol()
                W.symbols[i+1] = perm_autsymbol(p1*p2)
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
            reduced = simplify_perms!(W) && free_reduce!(W)
            deleteat!(W.symbols, find(x -> x.pow == 0, W.symbols))
        end
    end

    W.modified = false
    W.savedhash = hash(W.symbols,hash(typeof(W)))
    return W
end
