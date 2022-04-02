module MatrixGroups

using StaticArrays

using GroupsCore
using Groups
using KnuthBendix

import LinearAlgebra # Identity matrix
import Random # GroupsCore rand

export SpecialLinearGroup, SymplecticGroup

include("abstract.jl")

include("SLn.jl")
include("Spn.jl")

end # module
