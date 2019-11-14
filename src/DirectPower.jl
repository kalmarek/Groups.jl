export DirectPowerGroup, DirectPowerGroupElem

###############################################################################
#
#   DirectPowerGroup / DirectPowerGroupElem Constructors
#
###############################################################################

@doc doc"""
    DirectPowerGroup(G::Group, n::Int) <: Group
Implements `n`-fold direct product of `G`. The group operation is
`*` distributed component-wise, with component-wise identity as neutral element.
"""
struct DirectPowerGroup{N, T<:Group} <: Group
   group::T
end

DirectPowerGroup(G::Gr, N::Int) where Gr<:Group = DirectPowerGroup{N,Gr}(G)

function DirectPower(G::Group, H::Group)
   G == H || throw(DomainError(
      "Direct Powers are defined only for the same groups"))
   return DirectPowerGroup(G,2)
end

DirectPower(H::Group, G::DirectPowerGroup) = DirectPower(G,H)

function DirectPower(G::DirectPowerGroup{N}, H::Group) where N
   G.group == H || throw(DomainError(
      "Direct Powers are defined only for the same groups"))
   return DirectPowerGroup(G.group, N+1)
end

function DirectPower(R::AbstractAlgebra.Ring, n::Int)
   @warn "Creating DirectPower of the multilplicative group!"
   return DirectPowerGroup(MultiplicativeGroup(R), n)
end

struct DirectPowerGroupElem{N, T<:GroupElem} <: GroupElem
   elts::NTuple{N,T}
end

function DirectPowerGroupElem(v::Vector{GrEl}) where GrEl<:GroupElem
   return DirectPowerGroupElem(tuple(v...))
end

###############################################################################
#
#   Type and parent object methods
#
###############################################################################

elem_type(::Type{DirectPowerGroup{N,T}}) where {N,T} =
   DirectPowerGroupElem{N, elem_type(T)}

parent_type(::Type{DirectPowerGroupElem{N,T}}) where {N,T} =
   DirectPowerGroup{N, parent_type(T)}

parent(g::DirectPowerGroupElem{N, T}) where {N,T} =
   DirectPowerGroup(parent(first(g.elts)), N)

###############################################################################
#
#   AbstractVector interface
#
###############################################################################

size(g::DirectPowerGroupElem{N}) where N = (N,)
Base.IndexStyle(::Type{DirectPowerGroupElem}) = Base.LinearFast()
Base.getindex(g::DirectPowerGroupElem, i::Int) = g.elts[i]

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
function (G::DirectPowerGroup{N})(a::Vector, check::Bool=true) where N
   if check
      N == length(a) || throw(DomainError(
         "Can not coerce to DirectPowerGroup: lengths differ"))
      a = (G.group).(a)
   end
   return DirectPowerGroupElem(a)
end

function (G::DirectPowerGroup{N})(a::NTuple{N, GrEl}) where {N, GrEl}
   return DirectPowerGroupElem(G.group.(a))
end

(G::DirectPowerGroup{N})(a::Vararg{GrEl, N}) where {N, GrEl} = DirectPowerGroupElem(G.group.(a))

function Base.one(G::DirectPowerGroup{N}) where N
   return DirectPowerGroupElem(ntuple(i->one(G.group),N))
end

(G::DirectPowerGroup)(g::DirectPowerGroupElem) = G(g.elts)

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function hash(G::DirectPowerGroup{N}, h::UInt) where N
   return hash(G.group, hash(N, hash(DirectPowerGroup,h)))
end

function hash(g::DirectPowerGroupElem, h::UInt)
   return hash(g.elts, hash(DirectPowerGroupElem, h))
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, G::DirectPowerGroup{N}) where N
   print(io, "$(N)-fold direct product of $(G.group)")
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
function (==)(G::DirectPowerGroup{N}, H::DirectPowerGroup{M}) where {N,M}
   N == M || return false
   G.group == H.group || return false
   return true
end

@doc doc"""
    ==(g::DirectPowerGroupElem, h::DirectPowerGroupElem)
> Checks if two direct product group elements are the same.
"""
(==)(g::DirectPowerGroupElem, h::DirectPowerGroupElem) = g.elts == h.elts

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
function *(g::DirectPowerGroupElem{N}, h::DirectPowerGroupElem{N}, check::Bool=true) where N
   if check
      parent(g) == parent(h) || throw(DomainError(
         "Can not multiply elements of different groups!"))
   end
   return DirectPowerGroupElem(ntuple(i-> g.elts[i]*h.elts[i], N))
end

^(g::DirectPowerGroupElem, n::Integer) = Base.power_by_squaring(g, n)

@doc doc"""
    inv(g::DirectPowerGroupElem)
> Return the inverse of the given element in the direct product group.
"""
function inv(g::DirectPowerGroupElem{N}) where {N}
   return DirectPowerGroupElem(ntuple(i-> inv(g.elts[i]), N))
end

###############################################################################
#
#   Misc
#
###############################################################################

order(G::DirectPowerGroup{N}) where N = order(G.group)^N
length(G::DirectPowerGroup) = order(G)

function iterate(G::DirectPowerGroup{N}) where N
   elts = collect(G.group)

   indices = CartesianIndices(ntuple(i -> order(G.group), N))
   idx, s = iterate(indices)
   g = DirectPowerGroupElem(ntuple(i -> elts[idx[i]], N))
   return g, (elts, indices, s)
end

function iterate(G::DirectPowerGroup{N}, state) where N
   elts, indices, s = state
   res = iterate(indices, s)
   if res == nothing
      return nothing
   else
      idx, s = res
   end
   g = DirectPowerGroupElem(ntuple(i -> elts[idx[i]], N))
   return g, (elts, indices, s)
end

eltype(::Type{DirectPowerGroup{N, G}}) where {N, G} = DirectPowerGroupElem{N, elem_type(G)}
