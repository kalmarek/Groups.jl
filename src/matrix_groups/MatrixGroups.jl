module MatrixGroups

import LinearAlgebra # Identity matrix

using StaticArrays

using GroupsCore
import GroupsCore.Random # GroupsCore rand
using ..Groups
using Groups.KnuthBendix

export MatrixGroup, SpecialLinearGroup, SymplecticGroup

include("abstract.jl")

include("matrix_group.jl")
include("SLn.jl")
include("Spn.jl")

end # module
