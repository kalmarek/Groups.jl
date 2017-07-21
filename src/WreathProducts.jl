export WreathProduct, WreathProductElem

###############################################################################
#
#   WreathProduct / WreathProductElem
#
###############################################################################

doc"""
    WreathProduct <: Group
> Implements Wreath product of a group N by permutation (sub)group P < Sₖ,
> usually written as $N \wr P$.
> The multiplication inside wreath product is defined as
>    (n, σ) * (m, τ) = (n*ψ(σ)(m), σ*τ),
> where ψ:P → Aut(Nᵏ) is the permutation representation of Sₖ restricted to P.

# Arguments:
* `::Group` : the single factor of group N
* `::PermutationGroup` : full PermutationGroup
"""
immutable WreathProduct{T<:Group} <: Group
   N::DirectProductGroup{T}
   P::PermGroup

   function WreathProduct(G::Group, P::PermGroup)
      N = DirectProductGroup(G, P.n)
      return new(N, P)
   end
end

immutable WreathProductElem{T<:GroupElem} <: GroupElem
   n::DirectProductGroupElem{T}
   p::perm
   # parent::WreathProduct

   function WreathProductElem(n::DirectProductGroupElem, p::perm)
      length(n.elts) == parent(p).n || throw("Can't form WreathProductElem: lengths differ")
      return new(n, p)
   end
end

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

elem_type{T<:Group}(G::WreathProduct{T}) =
   WreathProductElem{G.N.n, elem_type(T)}

parent_type{T<:GroupElem}(::WreathProductElem{T}) =
   WreathProduct{parent_type(T)}

parent(g::WreathProductElem) = g.parent

###############################################################################
#
#   WreathProduct / WreathProductElem constructors
#
###############################################################################

# converts???

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (G::WreathProduct)(g::WreathProductElem)
   try
      G.N(g.n)
   catch
      throw("Can't coerce $(g.n) to $(G.N) factor of $G")
   end
   try
      G.P(g.p)
   catch
      throw("Can't coerce $(g.p) to $(G.P) factor of $G")
   end
   elt = WreathProductElem(G.N(g.n), G.P(g.p))
   elt.parent = G
   return elt
end

doc"""
    (G::WreathProduct)(n::DirectProductGroupElem, p::perm)
> Creates an element of wreath product `G` by coercing `n` and `p` to `G.N` and
> `G.P`, respectively.

"""
function (G::WreathProduct)(n::DirectProductGroupElem, p::perm)
   result = WreathProductElem(n,p)
   result.parent = G
   return result
end

(G::WreathProduct)() = G(G.N(), G.P())

doc"""
    (G::WreathProduct)(p::perm)
> Returns the image of permutation `p` in `G` via embedding `p -> (id,p)`.

"""
(G::WreathProduct)(p::perm) = G(G.N(), p)

doc"""
    (G::WreathProduct)(n::DirectProductGroupElem)
> Returns the image of `n` in `G` via embedding `n -> (n,())`. This is the
> embedding that makes sequence `1 -> N -> G -> P -> 1` exact.

"""
(G::WreathProduct)(n::DirectProductGroupElem) = G(n, G.P())

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function deepcopy_internal(g::WreathProductElem, dict::ObjectIdDict)
   G = parent(g)
   return G(deepcopy(g.n), deepcopy(g.p))
end

function hash(G::WreathProduct, h::UInt)
   return hash(G.N, hash(G.P, hash(WreathProduct, h)))
end

function hash(g::WreathProductElem, h::UInt)
   return hash(g.n, hash(g.p, hash(parent(g), h)))
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::WreathProduct)
   print(io, "Wreath Product of $(G.N.factors[1]) and $(G.P)")
end

function show(io::IO, g::WreathProductElem)
   # println(io, "Element of WreathProduct over $T of size $(size(X)):")
   # show(io, "text/plain", matrix_repr(X))
   print(io, "($(g.n)≀$(g.p))")
end

###############################################################################
#
#   Comparison
#
###############################################################################

function (==)(G::WreathProduct, H::WreathProduct)
   G.N == H.N || return false
   G.P == H.P || return false
   return true
end

function (==)(g::WreathProductElem, h::WreathProductElem)
   parent(g) == parent(h) || return false
   g.n == h.n || return false
   g.p == h.p || return false
   return true
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function wreath_multiplication(g::WreathProductElem, h::WreathProductElem)
   parent(g) == parent(h) || throw("Can not multiply elements from different
      groups!")
   G = parent(g)
   w = G.N((h.n).elts[inv(g.p).d])
   return G(g.n*w, g.p*h.p)
end

doc"""
    *(g::WreathProductElem, h::WreathProductElem)
> Return the wreath product group operation of elements, i.e.
>
> g*h = (g.n*g.p(h.n), g.p*h.p),
>
> where g.p(h.n) denotes the action of `g.p::perm` on
> `h.n::DirectProductGroupElem` via standard permutation of coordinates.
"""
(*)(g::WreathProductElem, h::WreathProductElem) = wreath_multiplication(g,h)


###############################################################################
#
#   Inversion
#
###############################################################################

doc"""
    inv(g::WreathProductElem)
> Returns the inverse of element of a wreath product, according to the formula
>   g^-1 = (g.n, g.p)^-1 = (g.p^-1(g.n^-1), g.p^-1).
"""
function inv(g::WreathProductElem)
   G = parent(g)
   w = G.N(inv(g.n).elts[g.p.d])
   return G(w, inv(g.p))
end

###############################################################################
#
#   Misc
#
###############################################################################

matrix_repr(g::WreathProductElem) = Any[matrix_repr(g.p) g.n]

function elements(G::WreathProduct)
   iter = Base.product(collect(elements(G.N)), collect(elements(G.P)))
   return (G(n)*G(p) for (n,p) in iter)
end

order(G::WreathProduct) = order(G.P)*order(G.N)
