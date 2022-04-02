module Groups

import Folds

using GroupsCore
import GroupsCore.Random

import OrderedCollections: OrderedSet

import KnuthBendix
import KnuthBendix: AbstractWord, Alphabet, Word
import KnuthBendix: alphabet

export MatrixGroups

export Alphabet, AutomorphismGroup, FreeGroup, FreeGroup, FPGroup, FPGroupElement, SpecialAutomorphismGroup, Homomorphism

export alphabet, evaluate, word, gens

# general constructions
include(joinpath("constructions", "constructions.jl"))
import .Constructions

include("types.jl")
include("hashing.jl")
include("normalform.jl")
include("autgroups.jl")
include("homomorphisms.jl")

include("aut_groups/sautFn.jl")
include("aut_groups/mcg.jl")

include("matrix_groups/MatrixGroups.jl")
using .MatrixGroups

include("abelianize.jl")

include("wl_ball.jl")
end # of module Groups
