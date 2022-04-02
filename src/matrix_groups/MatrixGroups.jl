module MatrixGroups

using GroupsCore
using Groups
using KnuthBendix

using LinearAlgebra # Identity matrix
using Random # GroupsCore rand

export SpecialLinearGroup, SymplecticGroup

include("abstract.jl")

include("SLn.jl")
include("Spn.jl")

end # module
