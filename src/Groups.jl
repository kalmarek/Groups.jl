module Groups

using GroupsCore
using ThreadsX
import KnuthBendix
import OrderedCollections: OrderedSet

export AutomorphismGroup, FreeGroup, FreeGroup, FPGroup, FPGroupElement, SpecialAutomorphismGroup
export alphabet, evaluate, word

include("new_types.jl")
include("new_hashing.jl")
include("normalform.jl")
include("new_autgroups.jl")

include("groups/sautFn.jl")

include("wl_ball.jl")
end # of module Groups
