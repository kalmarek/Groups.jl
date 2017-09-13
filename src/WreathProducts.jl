export WreathProduct, WreathProductElem

###############################################################################
#
#   WreathProduct / WreathProductElem
#
###############################################################################

doc"""
    WreathProduct{T<:Group} <: Group
> Implements Wreath product of a group $N$ by permutation (sub)group $P < S_k$,
> usually written as $N \wr P$.
> The multiplication inside wreath product is defined as
>    $$(n, \sigma) * (m, \tau) = (n\psi(\sigma)(m), \sigma\tau),$$
> where $\psi:P → Aut(N^k)$ is the permutation representation of $S_k$
> restricted to $P$.

# Arguments:
* `::Group` : the single factor of group $N$
* `::PermGroup` : full `PermutationGroup`
"""
immutable WreathProduct{T<:Group} <: Group
   N::DirectProductGroup{T}
   P::PermGroup

   function WreathProduct{T}(G::T, P::PermGroup) where {T}
      N = DirectProductGroup(G, P.n)
      return new(N, P)
   end
end

immutable WreathProductElem{T<:GroupElem} <: GroupElem
   n::DirectProductGroupElem{T}
   p::perm
   # parent::WreathProduct

   function WreathProductElem{T}(n::DirectProductGroupElem{T}, p::perm,
      check::Bool=true) where {T}
      if check
         length(n.elts) == parent(p).n || throw("Can't form WreathProductElem: lengths differ")
      end
      return new{T}(n, p)
   end
end

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

elem_type{T<:Group}(::WreathProduct{T}) = WreathProductElem{elem_type(T)}

parent_type{T<:GroupElem}(::Type{WreathProductElem{T}}) =
   WreathProduct{parent_type(T)}

parent(g::WreathProductElem) = WreathProduct(parent(g.n[1]), parent(g.p))

###############################################################################
#
#   WreathProduct / WreathProductElem constructors
#
###############################################################################

WreathProduct(G::Gr, P::PermGroup) where {Gr} = WreathProduct{Gr}(G, P)

WreathProductElem(n::DirectProductGroupElem{T}, p, check=true) where {T} = WreathProductElem{T}(n, p, check)

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (G::WreathProduct)(g::WreathProductElem)
   n = try
      G.N(g.n)
   catch
      throw("Can't coerce $(g.n) to $(G.N) factor of $G")
   end
   p = try
      G.P(g.p)
   catch
      throw("Can't coerce $(g.p) to $(G.P) factor of $G")
   end
   return WreathProductElem(n, p)
end

doc"""
    (G::WreathProduct)(n::DirectProductGroupElem, p::perm)
> Creates an element of wreath product `G` by coercing `n` and `p` to `G.N` and
> `G.P`, respectively.
"""
(G::WreathProduct)(n::DirectProductGroupElem, p::perm) = WreathProductElem(n,p)

(G::WreathProduct)() = WreathProductElem(G.N(), G.P(), false)

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
   return WreathProductElem(deepcopy(g.n), deepcopy(g.p), false)
end

function hash(G::WreathProduct, h::UInt)
   return hash(G.N, hash(G.P, hash(WreathProduct, h)))
end

function hash(g::WreathProductElem, h::UInt)
   return hash(g.n, hash(g.p, hash(WreathProductElem, h)))
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::WreathProduct)
   print(io, "Wreath Product of $(G.N.group) by $(G.P)")
end

function show(io::IO, g::WreathProductElem)
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
   g.n == h.n || return false
   g.p == h.p || return false
   return true
end

###############################################################################
#
#   Group operations
#
###############################################################################

doc"""
    *(g::WreathProductElem, h::WreathProductElem)
> Return the wreath product group operation of elements, i.e.
>
> `g*h = (g.n*g.p(h.n), g.p*h.p)`,
>
> where `g.p(h.n)` denotes the action of `g.p::perm` on
> `h.n::DirectProductGroupElem` via standard permutation of coordinates.
"""
function *(g::WreathProductElem, h::WreathProductElem)
   w = DirectProductGroupElem((h.n).elts[inv(g.p).d])
   return WreathProductElem(g.n*w, g.p*h.p, false)
end

doc"""
    inv(g::WreathProductElem)
> Returns the inverse of element of a wreath product, according to the formula
>   `g^-1 = (g.n, g.p)^-1 = (g.p^-1(g.n^-1), g.p^-1)`.
"""
function inv(g::WreathProductElem)
   w = DirectProductGroupElem(inv(g.n).elts[g.p.d])
   return WreathProductElem(w, inv(g.p), false)
end

###############################################################################
#
#   Misc
#
###############################################################################

matrix_repr(g::WreathProductElem) = Any[matrix_repr(g.p) g.n]

function elements(G::WreathProduct)
   iter = Base.product(collect(elements(G.N)), collect(elements(G.P)))
   return (WreathProductElem(n, p, false) for (n,p) in iter)
end

order(G::WreathProduct) = order(G.P)*order(G.N)
