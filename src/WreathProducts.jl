export WreathProduct, WreathProductElem

###############################################################################
#
#   WreathProduct / WreathProductElem
#
###############################################################################

doc"""
    WreathProduct(N, P) <: Group
> Implements Wreath product of a group `N` by permutation group $P = S_n$,
> usually written as $N \wr P$.
> The multiplication inside wreath product is defined as
>  >  `(n, σ) * (m, τ) = (n*σ(m), στ)`
> where `σ(m)` denotes the action (from the right) of the permutation group on
> `n-tuples` of elements from `N`

# Arguments:
* `N::Group` : the single factor of group $N$
* `P::Generic.PermGroup` : full `PermutationGroup`
"""
struct WreathProduct{T<:Group, I<:Integer} <: Group
   N::DirectProductGroup{T}
   P::Generic.PermGroup{I}

   function WreathProduct{T, I}(Gr::T, P::Generic.PermGroup{I}) where {T, I}
      N = DirectProductGroup(Gr, Int(P.n))
      return new(N, P)
   end
end

struct WreathProductElem{T<:GroupElem, I<:Integer} <: GroupElem
   n::DirectProductGroupElem{T}
   p::Generic.perm{I}
   # parent::WreathProduct

   function WreathProductElem{T, I}(n::DirectProductGroupElem{T}, p::Generic.perm{I},
      check::Bool=true) where {T, I}
      if check
         length(n.elts) == length(p) || throw(DomainError(
            "Can't form WreathProductElem: lengths differ"))
      end
      return new(n, p)
   end
end

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

elem_type(::Type{WreathProduct{T, I}}) where {T, I} = WreathProductElem{elem_type(T), I}

parent_type(::Type{WreathProductElem{T, I}}) where {T, I} =
   WreathProduct{parent_type(T), I}

parent(g::WreathProductElem) = WreathProduct(parent(g.n[1]), parent(g.p))

###############################################################################
#
#   WreathProduct / WreathProductElem constructors
#
###############################################################################

WreathProduct(G::T, P::Generic.PermGroup{I}) where {T, I} = WreathProduct{T, I}(G, P)

WreathProduct(G::T, P::Generic.PermGroup{I}) where {T<:AbstractAlgebra.Ring, I} = WreathProduct(AddGrp(G), P)

WreathProductElem(n::DirectProductGroupElem{T}, p::Generic.perm{I}, check=true) where {T,I} = WreathProductElem{T,I}(n, p, check)

WreathProductElem(n::DirectProductGroupElem{T}, p::Generic.perm{I}, check=true) where {T<:AbstractAlgebra.RingElem, I} = WreathProductElem(DirectProductGroupElem(AddGrpElem.(n.elts)), p, check)

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (G::WreathProduct)(g::WreathProductElem)
   n = try
      G.N(g.n)
   catch
      throw(DomainError("Can't coerce $(g.n) to $(G.N) factor of $G"))
   end
   p = try
      G.P(g.p)
   catch
      throw(DomainError("Can't coerce $(g.p) to $(G.P) factor of $G"))
   end
   return WreathProductElem(n, p)
end

doc"""
    (G::WreathProduct)(n::DirectProductGroupElem, p::Generic.perm)
> Creates an element of wreath product `G` by coercing `n` and `p` to `G.N` and
> `G.P`, respectively.
"""
(G::WreathProduct)(n::DirectProductGroupElem, p::Generic.perm) = WreathProductElem(n,p)

(G::WreathProduct)() = WreathProductElem(G.N(), G.P(), false)

doc"""
    (G::WreathProduct)(p::Generic.perm)
> Returns the image of permutation `p` in `G` via embedding `p -> (id,p)`.
"""
(G::WreathProduct)(p::Generic.perm) = G(G.N(), p)

doc"""
    (G::WreathProduct)(n::DirectProductGroupElem)
> Returns the image of `n` in `G` via embedding `n -> (n,())`. This is the
> embedding that makes sequence `1 -> N -> G -> P -> 1` exact.
"""
(G::WreathProduct)(n::DirectProductGroupElem) = G(n, G.P())

(G::WreathProduct)(n,p) = G(G.N(n), G.P(p))

###############################################################################
#
#   Basic manipulation
#
###############################################################################

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

(p::perm)(n::DirectProductGroupElem) = DirectProductGroupElem(n.elts[p.d])

doc"""
    *(g::WreathProductElem, h::WreathProductElem)
> Return the wreath product group operation of elements, i.e.
>
> `g*h = (g.n*g.p(h.n), g.p*h.p)`,
>
> where `g.p(h.n)` denotes the action of `g.p::Generic.perm` on
> `h.n::DirectProductGroupElem` via standard permutation of coordinates.
"""
function *(g::WreathProductElem, h::WreathProductElem)
   return WreathProductElem(g.n*g.p(h.n), g.p*h.p, false)
end

doc"""
    inv(g::WreathProductElem)
> Returns the inverse of element of a wreath product, according to the formula
>   `g^-1 = (g.n, g.p)^-1 = (g.p^-1(g.n^-1), g.p^-1)`.
"""
function inv(g::WreathProductElem)
   pinv = inv(g.p)
   return WreathProductElem(pinv(inv(g.n)), pinv, false)
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
