import Base: ×

export DirectProductGroup, DirectProductGroupElem
export MultiplicativeGroup, MltGrp, MltGrpElem
export AdditiveGroup, AddGrp, AddGrpElem

###############################################################################
#
#   MltGrp/MltGrpElem & AddGrp/AddGrpElem
#   (a thin wrapper for multiplicative/additive group of a Ring)
#
###############################################################################

for (Gr, Elem) in [(:MltGrp, :MltGrpElem), (:AddGrp, :AddGrpElem)]
   @eval begin
      struct $Gr{T<:AbstractAlgebra.Ring} <: AbstractAlgebra.Group
         obj::T
      end

      struct $Elem{T<:AbstractAlgebra.RingElem} <: AbstractAlgebra.GroupElem
         elt::T
      end

      ==(g::$Elem, h::$Elem) = g.elt == h.elt
      ==(G::$Gr, H::$Gr) = G.obj == H.obj

      elem_type(::Type{$Gr{T}}) where T = $Elem{elem_type(T)}
      parent_type(::Type{$Elem{T}}) where T = $Gr{parent_type(T)}
      parent(g::$Elem) = $Gr(parent(g.elt))
   end
end

MultiplicativeGroup = MltGrp
AdditiveGroup = AddGrp

(G::MltGrp)(g::MltGrpElem) = MltGrpElem(G.obj(g.elt))

function (G::MltGrp)(g)
   r = (G.obj)(g)
   isunit(r) || throw(ArgumentError("Cannot coerce to multplicative group: $r is not invertible!"))
   return MltGrpElem(r)
end

(G::AddGrp)(g) = AddGrpElem((G.obj)(g))

(G::MltGrp)() = MltGrpElem(G.obj(1))
(G::AddGrp)() = AddGrpElem(G.obj())

inv(g::MltGrpElem) = MltGrpElem(inv(g.elt))
inv(g::AddGrpElem) = AddGrpElem(-g.elt)

for (Elem, op) in ([:MltGrpElem, :*], [:AddGrpElem, :+])
   @eval begin

      ^(g::$Elem, n::Integer) = $Elem(op(g.elt, n))

      function *(g::$Elem, h::$Elem)
         parent(g) == parent(h) || throw("Cannot multiply elements of different parents")
         return $Elem($op(g.elt,h.elt))
      end
   end
end

Base.show(io::IO, G::MltGrp) = print(io, "The multiplicative group of $(G.obj)")
Base.show(io::IO, G::AddGrp) = print(io, "The additive group of $(G.obj)")

Base.show(io::IO, g::Union{MltGrpElem, AddGrpElem}) = show(io, g.elt)

gens(F::AbstractAlgebra.Field) = elem_type(F)[gen(F)]

order(G::AddGrp{<:AbstractAlgebra.GFField}) = order(G.obj)
elements(G::AddGrp{F}) where F <: AbstractAlgebra.GFField = (G((i-1)*G.obj(1)) for i in 1:order(G))

order(G::MltGrp{<:AbstractAlgebra.GFField}) = order(G.obj) - 1
elements(G::MltGrp{F}) where F <: AbstractAlgebra.GFField = (G(i*G.obj(1)) for i in 1:order(G))


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

struct DirectProductGroup{T<:Group} <: Group
   group::T
   n::Int
end

struct DirectProductGroupElem{T<:GroupElem} <: GroupElem
   elts::Vector{T}
end

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

elem_type(G::DirectProductGroup{T}) where {T} =
   DirectProductGroupElem{elem_type(G.group)}

parent_type(::Type{DirectProductGroupElem{T}}) where {T} =
   DirectProductGroup{parent_type(T)}

parent(g::DirectProductGroupElem) =
   DirectProductGroup(parent(first(g.elts)), length(g.elts))

###############################################################################
#
#   AbstractVector interface
#
###############################################################################

Base.size(g::DirectProductGroupElem) = size(g.elts)
Base.IndexStyle(::Type{DirectProductGroupElem}) = Base.LinearFast()
Base.getindex(g::DirectProductGroupElem, i::Int) = g.elts[i]
function Base.setindex!(g::DirectProductGroupElem{T}, v::T, i::Int) where {T}
   parent(v) == parent(first(g.elts)) || throw("$g is not an element of $i-th factor of $(parent(G))")
   g.elts[i] = v
   return g
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

DirectProductGroup(R::T, n::Int) where {T<:AbstractAlgebra.Ring} =
DirectProductGroup(AdditiveGroup(R), n)

function ×(G::DirectProductGroup{T}, H::Group) where T <: Union{AdditiveGroup, MultiplicativeGroup}
   G.group == T(H) || throw(ArgumentError("Direct products are defined only for the same groups"))
   return DirectProductGroup(G.group,G.n+1)
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

doc"""
    (G::DirectProductGroup)(a::Vector, check::Bool=true)
> Constructs element of the $n$-fold direct product group `G` by coercing each
> element of vector `a` to `G.group`. If `check` flag is set to `false` neither
> check on the correctness nor coercion is performed.
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

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::DirectProductGroup)
   print(io, "$(G.n)-fold direct product of $(G.group)")
end

function show(io::IO, g::DirectProductGroupElem)
   print(io, "[$(join(g.elts,","))]")
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
   g.elts == h.elts || return false
   return true
end

###############################################################################
#
#   Group operations
#
###############################################################################

doc"""
    *(g::DirectProductGroupElem, h::DirectProductGroupElem)
> Return the direct-product group operation of elements, i.e. component-wise
> operation as defined by `operations` field of the parent object.
"""
function *(g::DirectProductGroupElem{T}, h::DirectProductGroupElem{T}, check::Bool=true) where {T}
   if check
      parent(g) == parent(h) || throw("Can not multiply elements of different groups!")
   end
   return DirectProductGroupElem([a*b for (a,b) in zip(g.elts,h.elts)])
end

doc"""
    inv(g::DirectProductGroupElem)
> Return the inverse of the given element in the direct product group.
"""
function inv(g::DirectProductGroupElem{T}) where {T<:GroupElem}
   return DirectProductGroupElem([inv(a) for a in g.elts])
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
