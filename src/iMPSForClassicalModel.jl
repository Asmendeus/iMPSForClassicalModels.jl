module iMPSForClassicalModel

using Reexport
using TensorKit, OptimKit, KrylovKit
@reexport import Base: +, -, *, /, ==, promote_rule, convert, length, size, show, getindex, setindex!, lastindex, keys, similar, merge, merge!, iterate, complex

# MPS
export AbstractInfiniteMPS
export InfiniteMPS, iMPS
include("iMPS/AbstractInfiniteMPS.jl")
include("iMPS/InfiniteMPS.jl")

# MPO
export AbstractInfiniteMPO
export InfiniteMPO, iMPO
export MultirowInfiniteMPO, miMPO
include("iMPO/AbstractInfiniteMPO.jl")
include("iMPO/InfiniteMPO.jl")
include("iMPO/MultirowInfiniteMPO.jl")

# Environment

# Algorithm
include("Algorithm/iTEBD.jl")
include("Algorithm/VUMPS.jl")

end # module iMPSForClassicalModel
