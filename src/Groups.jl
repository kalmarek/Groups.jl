module Groups

using GroupsCore
using Folds
import KnuthBendix
import KnuthBendix: AbstractWord, Alphabet, Word
import KnuthBendix: alphabet
import Random

import OrderedCollections: OrderedSet

export Alphabet, AutomorphismGroup, FreeGroup, FreeGroup, FPGroup, FPGroupElement, SpecialAutomorphismGroup
export alphabet, evaluate, word, gens

include("types.jl")
include("hashing.jl")
include("normalform.jl")
include("autgroups.jl")

include("aut_groups/sautFn.jl")
include("aut_groups/mcg.jl")


include("wl_ball.jl")
end # of module Groups
