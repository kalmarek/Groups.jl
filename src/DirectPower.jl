export DirectPowerGroup, DirectPowerGroupElem
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
      eltype(::Type{$Gr{T}}) where T = $Elem{elem_type(T)}
      parent_type(::Type{$Elem{T}}) where T = $Gr{parent_type(T)}
      parent(g::$Elem) = $Gr(parent(g.elt))
   end
end

MultiplicativeGroup = MltGrp
AdditiveGroup = AddGrp

(G::MltGrp)(g::MltGrpElem) = MltGrpElem(G.obj(g.elt))

function (G::MltGrp)(g)
   r = (G.obj)(g)
   isunit(r) || throw(DomainError(
      "Cannot coerce to multplicative group: $r is not invertible!"))
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
         parent(g) == parent(h) || throw(DomainError(
            "Cannot multiply elements of different parents"))
         return $Elem($op(g.elt,h.elt))
      end
   end
end

show(io::IO, G::MltGrp) = print(io, "The multiplicative group of $(G.obj)")
show(io::IO, G::AddGrp) = print(io, "The additive group of $(G.obj)")

show(io::IO, g::Union{MltGrpElem, AddGrpElem}) = show(io, g.elt)

gens(F::AbstractAlgebra.Field) = elem_type(F)[gen(F)]

order(G::AddGrp{<:AbstractAlgebra.GFField}) = order(G.obj)
elements(G::AddGrp{F}) where F <: AbstractAlgebra.GFField = (G((i-1)*G.obj(1)) for i in 1:order(G))

order(G::MltGrp{<:AbstractAlgebra.GFField}) = order(G.obj) - 1
elements(G::MltGrp{F}) where F <: AbstractAlgebra.GFField = (G(i*G.obj(1)) for i in 1:order(G))

function iterate(G::AddGrp, s=0)
   if s >= order(G)
      return nothing
   else
      g, s = iterate(G.obj,s)
      return G(g), s
   end
end

function iterate(G::MltGrp, s=0)
   if s > order(G)
      return nothing
   else
      g, s = iterate(G.obj, s)
      if g == G.obj()
         g, s = iterate(G.obj, s)
      end
      return G(g), s
   end
end

###############################################################################
#
#   DirectPowerGroup / DirectPowerGroupElem
#
###############################################################################

@doc doc"""
    DirectPowerGroup(G::Group, n::Int) <: Group
Implements `n`-fold direct product of `G`. The group operation is
`*` distributed component-wise, with component-wise identity as neutral element.
"""
struct DirectPowerGroup{T<:Group} <: Group
   group::T
   n::Int
end

struct DirectPowerGroupElem{T<:GroupElem} <: GroupElem
   elts::Vector{T}
end

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

elem_type(::Type{DirectPowerGroup{T}}) where {T} =
   DirectPowerGroupElem{elem_type(T)}

parent_type(::Type{DirectPowerGroupElem{T}}) where {T} =
   DirectPowerGroup{parent_type(T)}

parent(g::DirectPowerGroupElem) =
   DirectPowerGroup(parent(first(g.elts)), length(g.elts))

###############################################################################
#
#   AbstractVector interface
#
###############################################################################

size(g::DirectPowerGroupElem) = size(g.elts)
Base.IndexStyle(::Type{DirectPowerGroupElem}) = Base.LinearFast()
Base.getindex(g::DirectPowerGroupElem, i::Int) = g.elts[i]

function Base.setindex!(g::DirectPowerGroupElem{T}, v::T, i::Int) where {T}
   parent(v) == parent(g.elts[i]) || throw(DomainError(
      "$g is not an element of $i-th factor of $(parent(G))"))
   g.elts[i] = v
   return g
end

function Base.setindex!(g::DirectPowerGroupElem{T}, v::S, i::Int) where {T, S}
   g.elts[i] = parent(g.elts[i])(v)
   return g
end

###############################################################################
#
#   DirectPowerGroup / DirectPowerGroupElem constructors
#
###############################################################################

function DirectPower(G::Group, H::Group)
   G == H || throw(DomainError(
      "Direct Powers are defined only for the same groups"))
   return DirectPowerGroup(G,2)
end

DirectPower(H::Group, G::DirectPowerGroup) = DirectPower(G,H)

function DirectPower(G::DirectPowerGroup, H::Group)
   G.group == H || throw(DomainError(
      "Direct products are defined only for the same groups"))
   return DirectPowerGroup(G.group,G.n+1)
end

function DirectPower(R::AbstractAlgebra.Ring, n::Int)
   @warn "Creating DirectPower of the multilplicative group!"
   return DirectPowerGroup(R, n)
end

###############################################################################
#
#   Parent object call overloads
#
###############################################################################

@doc doc"""
    (G::DirectPowerGroup)(a::Vector, check::Bool=true)
> Constructs element of the $n$-fold direct product group `G` by coercing each
> element of vector `a` to `G.group`. If `check` flag is set to `false` neither
> check on the correctness nor coercion is performed.
"""
function (G::DirectPowerGroup)(a::Vector, check::Bool=true)
   if check
      G.n == length(a) || throw(DomainError(
         "Can not coerce to DirectPowerGroup: lengths differ"))
      a = (G.group).(a)
   end
   return DirectPowerGroupElem(a)
end

(G::DirectPowerGroup)() = DirectPowerGroupElem([G.group() for _ in 1:G.n])

(G::DirectPowerGroup)(g::DirectPowerGroupElem) = G(g.elts)

(G::DirectPowerGroup)(a::Vararg{T, N}) where {T, N} = G([a...])

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function hash(G::DirectPowerGroup, h::UInt)
   return hash(G.group, hash(G.n, hash(DirectPowerGroup,h)))
end

function hash(g::DirectPowerGroupElem, h::UInt)
   return hash(g.elts, hash(parent(g), hash(DirectPowerGroupElem, h)))
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::DirectPowerGroup)
   print(io, "$(G.n)-fold direct product of $(G.group)")
end

function show(io::IO, g::DirectPowerGroupElem)
   print(io, "[$(join(g.elts,","))]")
end

###############################################################################
#
#   Comparison
#
###############################################################################

@doc doc"""
    ==(g::DirectPowerGroup, h::DirectPowerGroup)
> Checks if two direct product groups are the same.
"""
function (==)(G::DirectPowerGroup, H::DirectPowerGroup)
   G.group == H.group || return false
   G.n == G.n || return false
   return true
end

@doc doc"""
    ==(g::DirectPowerGroupElem, h::DirectPowerGroupElem)
> Checks if two direct product group elements are the same.
"""
function (==)(g::DirectPowerGroupElem, h::DirectPowerGroupElem)
   g.elts == h.elts || return false
   return true
end

###############################################################################
#
#   Group operations
#
###############################################################################

@doc doc"""
    *(g::DirectPowerGroupElem, h::DirectPowerGroupElem)
> Return the direct-product group operation of elements, i.e. component-wise
> operation as defined by `operations` field of the parent object.
"""
function *(g::DirectPowerGroupElem, h::DirectPowerGroupElem, check::Bool=true)
   if check
      parent(g) == parent(h) || throw(DomainError(
         "Can not multiply elements of different groups!"))
   end
   return DirectPowerGroupElem([a*b for (a,b) in zip(g.elts,h.elts)])
end

^(g::DirectPowerGroupElem, n::Integer) = Base.power_by_squaring(g, n)

@doc doc"""
    inv(g::DirectPowerGroupElem)
> Return the inverse of the given element in the direct product group.
"""
function inv(g::DirectPowerGroupElem{T}) where {T<:GroupElem}
   return DirectPowerGroupElem([inv(a) for a in g.elts])
end

###############################################################################
#
#   Misc
#
###############################################################################

struct DirectPowerIter{GrEl<:AbstractAlgebra.GroupElem}
   N::Int
   elts::Vector{GrEl}
   totalorder::Int
   orderG::Int
end

function DirectPowerIter(G::Gr, N::Integer) where {Gr<:AbstractAlgebra.Group}
   return DirectPowerIter{elem_type(G)}(N, collect(G), order(G)^N, order(G))
end

length(DPIter::DirectPowerIter) = DPIter.totalorder

function iterate(DPIter::DirectPowerIter, state=0)
   if state >= DPIter.totalorder
      return nothing
   end
   idx = Tuple(CartesianIndices(ntuple(i -> DPIter.orderG, DPIter.N))[state+1])
   return DirectPowerGroupElem([DPIter.elts[i] for i in idx]), state+1
end

eltype(::Type{DirectPowerIter{GrEl}}) where {GrEl} = DirectPowerGroupElem{GrEl}

@doc doc"""
    elements(G::DirectPowerGroup)
> Returns `generator` that produces all elements of group `G` (provided that
> `G.group` implements the `elements` method).
"""
elements(G::DirectPowerGroup) = DirectPowerIter(G.group, G.n)

@doc doc"""
    order(G::DirectPowerGroup)
> Returns the order (number of elements) in the group.
"""
order(G::DirectPowerGroup) = order(G.group)^G.n
