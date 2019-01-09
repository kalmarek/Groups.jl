export WreathProduct, WreathProductElem

###############################################################################
#
#   WreathProduct / WreathProductElem
#
###############################################################################

@doc doc"""
    WreathProduct(N, P) <: Group
> Implements Wreath product of a group `N` by permutation group $P = S_n$,
> usually written as $N \wr P$.
> The multiplication inside wreath product is defined as
>  >  `(n, σ) * (m, τ) = (n*σ(m), στ)`
> where `σ(m)` denotes the action (from the right) of the permutation group on
> `n-tuples` of elements from `N`

# Arguments:
* `N::Group` : the single factor of the group $N$
* `P::Generic.PermGroup` : full `PermutationGroup`
"""
struct WreathProduct{N, T<:Group, PG<:Generic.PermGroup} <: Group
   N::DirectPowerGroup{N, T}
   P::PG

   function WreathProduct(Gr::T, P::PG) where {T, PG<:Generic.PermGroup}
      N = DirectPowerGroup(Gr, Int(P.n))
      return new{Int(P.n), T, PG}(N, P)
   end
end

struct WreathProductElem{N, T<:GroupElem, P<:Generic.perm} <: GroupElem
   n::DirectPowerGroupElem{N, T}
   p::P

   function WreathProductElem(n::DirectPowerGroupElem{N,T}, p::P,
      check::Bool=true) where {N, T, P<:Generic.perm}
      if check
         N == length(p.d) || throw(DomainError(
            "Can't form WreathProductElem: lengths differ"))
      end
      return new{N, T, P}(n, p)
   end
end

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

elem_type(::Type{WreathProduct{N, T, PG}}) where {N, T, PG} = WreathProductElem{N, elem_type(T), elem_type(PG)}

parent_type(::Type{WreathProductElem{N, T, P}}) where {N, T, P} =
   WreathProduct{N, parent_type(T), parent_type(P)}

parent(g::WreathProductElem) = WreathProduct(parent(g.n[1]), parent(g.p))

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

function (G::WreathProduct{N})(g::WreathProductElem{N}) where {N}
   n = G.N(g.n)
   p = G.P(g.p)
   return WreathProductElem(n, p)
end

@doc doc"""
    (G::WreathProduct)(n::DirectPowerGroupElem, p::Generic.perm)
> Creates an element of wreath product `G` by coercing `n` and `p` to `G.N` and
> `G.P`, respectively.
"""
(G::WreathProduct)(n::DirectPowerGroupElem, p::Generic.perm) = WreathProductElem(n,p)

(G::WreathProduct)() = WreathProductElem(G.N(), G.P(), false)

@doc doc"""
    (G::WreathProduct)(p::Generic.perm)
> Returns the image of permutation `p` in `G` via embedding `p -> (id,p)`.
"""
(G::WreathProduct)(p::Generic.perm) = G(G.N(), p)

@doc doc"""
    (G::WreathProduct)(n::DirectPowerGroupElem)
> Returns the image of `n` in `G` via embedding `n -> (n,())`. This is the
> embedding that makes the sequence `1 -> N -> G -> P -> 1` exact.
"""
(G::WreathProduct)(n::DirectPowerGroupElem) = G(n, G.P())

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

(p::perm)(n::DirectPowerGroupElem) = DirectPowerGroupElem(n.elts[p.d])

@doc doc"""
    *(g::WreathProductElem, h::WreathProductElem)
> Return the wreath product group operation of elements, i.e.
>
> `g*h = (g.n*g.p(h.n), g.p*h.p)`,
>
> where `g.p(h.n)` denotes the action of `g.p::Generic.perm` on
> `h.n::DirectPowerGroupElem` via standard permutation of coordinates.
"""
function *(g::WreathProductElem, h::WreathProductElem)
   return WreathProductElem(g.n*g.p(h.n), g.p*h.p, false)
end

^(g::WreathProductElem, n::Integer) = Base.power_by_squaring(g, n)

@doc doc"""
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

function iterate(G::WreathProduct)
   n, state_N = iterate(G.N)
   p, state_P = iterate(G.P)
   return G(n,p), (state_N, p, state_P)
end

function iterate(G::WreathProduct, state)
   state_N, p, state_P = state
   res = iterate(G.N, state_N)
   
   if res == nothing
      resP = iterate(G.P, state_P)
      if resP == nothing
         return nothing
      else
         n, state_N = iterate(G.N)
         p, state_P = resP
      end
   else
      n, state_N = res
   end
   
   return G(n,p), (state_N, p, state_P)
end

eltype(::Type{WreathProduct{N,G,PG}}) where {N,G,PG} = WreathProductElem{N, elem_type(G), elem_type(PG)}

order(G::WreathProduct) = order(G.P)*order(G.N)
length(G::WreathProduct) = order(G)
