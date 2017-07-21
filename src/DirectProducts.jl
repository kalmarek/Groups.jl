import Base: ×

export DirectProductGroup, DirectProductGroupElem

###############################################################################
#
#   DirectProductGroup / DirectProductGroupElem
#
###############################################################################

doc"""
   DirectProductGroup(G::Group, n::Int) <: Group
Implements `n`-fold direct product of `G`. The group operation is
`*` distributed component-wise, with component-wise identity as neutral element.
"""

immutable DirectProductGroup{T<:Group} <: Group
   group::T
   n::Int
end

immutable DirectProductGroupElem{T<:GroupElem} <: GroupElem
   elts::Vector{T}
end

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

elem_type{T<:Group}(G::DirectProductGroup{T}) =
   DirectProductGroupElem{elem_type(G.group)}

parent_type{T<:GroupElem}(::Type{DirectProductGroupElem{T}}) =
   DirectProductGroup{parent_type(T)}

parent(g::DirectProductGroupElem) =
   DirectProductGroup(parent(first(g.elts)), length(g.elts))

###############################################################################
#
#   AbstractVector interface
#
###############################################################################

Base.size(g::DirectProductGroupElem) = size(g.elts)
Base.linearindexing(::Type{DirectProductGroupElem}) = Base.LinearFast()
Base.getindex(g::DirectProductGroupElem, i::Int) = g.elts[i]
function Base.setindex!{T<:GroupElem}(g::DirectProductGroupElem{T}, v::T, i::Int)
   p.part[i] = v
   return p
end

###############################################################################
#
#   DirectProductGroup / DirectProductGroupElem constructors
#
###############################################################################

function ×(G::Group, H::Group)
   G == H || throw("Direct products are defined only for the same groups")
   return DirectProductGroup(G,2)
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

doc"""
    (G::DirectProductGroup)(a::Vector, check::Bool=true)
> Constructs element of the $n$-fold direct product group `G` by coercing each
> element of vector `a` to `G.group`. If `check` flag is set to `false` no
> checks on the correctness are performed.
"""
function (G::DirectProductGroup)(a::Vector, check::Bool=true)
   if check
      G.n == length(a) || throw("Can not coerce to DirectProductGroup: lengths differ")
      a = G.group.(a)
   end
   return DirectProductGroupElem(a)
end

(G::DirectProductGroup)() = DirectProductGroupElem([G.group() for _ in 1:G.n])

(G::DirectProductGroup)(g::DirectProductGroupElem) = G(g.elts)

(G::DirectProductGroup){T<:GroupElem, N}(a::Vararg{T, N}) = G([a...])

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function deepcopy_internal(g::DirectProductGroupElem, dict::ObjectIdDict)
   G = parent(g)
   return G(deepcopy(g.elts))
end

function hash(G::DirectProductGroup, h::UInt)
   return hash(G.group, hash(G.n, hash(DirectProductGroup,h)))
end

function hash(g::DirectProductGroupElem, h::UInt)
   return hash(g.elts, hash(parent(g), hash(DirectProductGroupElem, h)))
end

doc"""
    eye(G::DirectProductGroup)
> Return the identity element for the given direct product of groups.

"""
eye(G::DirectProductGroup) = G()

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::DirectProductGroup)
    println(io, "$(G.n)-fold direct product of $(G.group)")
end

function show(io::IO, g::DirectProductGroupElem)
    print(io, "($(join(g.elts,",")))")
end

###############################################################################
#
#   Comparison
#
###############################################################################

doc"""
    ==(g::DirectProductGroup, h::DirectProductGroup)
> Checks if two direct product groups are the same.
"""
function (==)(G::DirectProductGroup, H::DirectProductGroup)
    G.group == H.group || return false
    G.n == G.n || return false
    return true
end

doc"""
    ==(g::DirectProductGroupElem, h::DirectProductGroupElem)
> Checks if two direct product group elements are the same.
"""
function (==)(g::DirectProductGroupElem, h::DirectProductGroupElem)
    parent(g) == parent(h) || return false
    g.elts == h.elts || return false
    return true
end

###############################################################################
#
#   Binary operators
#
###############################################################################

doc"""
    *(g::DirectProductGroupElem, h::DirectProductGroupElem)
> Return the direct-product group operation of elements, i.e. component-wise
> operation as defined by `operations` field of the parent object.
"""
# TODO: dirty hack around `+/*` operations
function *{T<:GroupElem}(g::DirectProductGroupElem{T}, h::DirectProductGroupElem{T}, check::Bool=true)
   if check
      parent(g) == parent(h) || throw("Can not multiply elements of different groups!")
   end
   return DirectProductGroupElem([a*b for (a,b) in zip(g.elts,h.elts)])
end

function *{T<:RingElem}(g::DirectProductGroupElem{T}, h::DirectProductGroupElem{T}, check::Bool=true)
   if check
      parent(g) == parent(h) || throw("Can not multiply elements of different groups!")
   end
   return DirectProductGroupElem([a+b for (a,b) in zip(g.elts,h.elts)])
end

###############################################################################
#
#   Inversion
#
###############################################################################

doc"""
    inv(g::DirectProductGroupElem)
> Return the inverse of the given element in the direct product group.
"""
# TODO: dirty hack around `+/*` operation
function inv{T<:GroupElem}(g::DirectProductGroupElem{T})
   return DirectProductGroupElem([inv(a) for a in g.elts])
end

function inv{T<:RingElem}(g::DirectProductGroupElem{T})
   return DirectProductGroupElem([-a for a in g.elts])
end

###############################################################################
#
#   Misc
#
###############################################################################

doc"""
    elements(G::DirectProductGroup)
> Returns `generator` that produces all elements of group `G` (provided that
> `G.group` implements the `elements` method).
"""
# TODO: can Base.product handle generators?
#   now it returns nothing's so we have to collect ellements...
function elements(G::DirectProductGroup)
    elts = collect(elements(G.group))
    cartesian_prod = Base.product([elts for _ in 1:G.n]...)
    return (DirectProductGroupElem([elt...]) for elt in cartesian_prod)
end

doc"""
    order(G::DirectProductGroup)
> Returns the order (number of elements) in the group.

"""
order(G::DirectProductGroup) = order(G.group)^G.n
