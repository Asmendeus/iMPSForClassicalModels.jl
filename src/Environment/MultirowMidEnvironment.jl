"""
    mutable struct MultirowMidEnvironment{N} <: AbstractEnvironment{N}
        const BC::AbstractVector{<:AdjointMPSTensor}
        const O::AbstractMatrix{MPOTensor}
        const AC::AbstractVector{<:MPSTensor}
    end

Wrapper type for Mid environment tensor's generating environment.

Graphic presentation:

    —— BC ——
       |
    —— OW ——    (W = N-2)
       |
       ⋮
       |
    —— O1 ——
       |
    —— AC ——

# Constructors
    MultirowMidEnvironment{N}(BC::AbstractVector{<:AdjointMPSTensor}, O::AbstractMatrix{MPOTensor}, AC::AbstractVector{<:MPSTensor})
    MultirowMidEnvironment(BC::AbstractVector{<:AdjointMPSTensor}, O::AbstractMatrix{MPOTensor}, AC::AbstractVector{<:MPSTensor})
    MultirowMidEnvironment{N}(AC::AbstractVector{<:MPSTensor}, O::AbstractMatrix{MPOTensor}, BC::AbstractVector{<:AdjointMPSTensor})
    MultirowMidEnvironment(AC::AbstractVector{<:MPSTensor}, O::AbstractMatrix{MPOTensor}, BC::AbstractVector{<:AdjointMPSTensor})
    MultirowMidEnvironment{N}(BC::AbstractVector{<:AdjointMPSTensor}, AC::AbstractVector{<:MPSTensor}, O::AbstractMatrix{MPOTensor})
    MultirowMidEnvironment(BC::AbstractVector{<:AdjointMPSTensor}, AC::AbstractVector{<:MPSTensor}, O::AbstractMatrix{MPOTensor})
    MultirowMidEnvironment{N}(AC::AbstractVector{<:MPSTensor}, BC::AbstractVector{<:AdjointMPSTensor}, O::AbstractMatrix{MPOTensor})
    MultirowMidEnvironment(AC::AbstractVector{<:MPSTensor}, BC::AbstractVector{<:AdjointMPSTensor}, O::AbstractMatrix{MPOTensor})
    MultirowMidEnvironment{N}(O::AbstractMatrix{MPOTensor}, BC::AbstractVector{<:AdjointMPSTensor}, AC::AbstractVector{<:MPSTensor})
    MultirowMidEnvironment(O::AbstractMatrix{MPOTensor}, BC::AbstractVector{<:AdjointMPSTensor}, AC::AbstractVector{<:MPSTensor})
    MultirowMidEnvironment{N}(O::AbstractMatrix{MPOTensor}, AC::AbstractVector{<:MPSTensor}, BC::AbstractVector{<:AdjointMPSTensor})
    MultirowMidEnvironment(O::AbstractMatrix{MPOTensor}, AC::AbstractVector{<:MPSTensor}, BC::AbstractVector{<:AdjointMPSTensor})
"""
mutable struct MultirowMidEnvironment{N} <: AbstractEnvironment{N}
    const BC::AbstractVector{<:AdjointMPSTensor}
    const O::AbstractMatrix{MPOTensor}
    const AC::AbstractVector{<:MPSTensor}

    function MultirowMidEnvironment{N}(BC::AbstractVector{<:AdjointMPSTensor}, O::AbstractMatrix{MPOTensor}, AC::AbstractVector{<:MPSTensor}) where N
        length(O) == N - 2 || throw(ArgumentError("The width of multirow MPO doesn't match width of `MultirowMidEnvironment`"))
        return new{N}(BC, O, AC)
    end
    function MultirowMidEnvironment(BC::AbstractVector{<:AdjointMPSTensor}, O::AbstractMatrix{MPOTensor}, AC::AbstractVector{<:MPSTensor})
        N = length(O) + 2
        return new{N}(BC, O, AC)
    end

    function MultirowMidEnvironment{N}(AC::AbstractVector{<:MPSTensor}, O::AbstractMatrix{MPOTensor}, BC::AbstractVector{<:AdjointMPSTensor}) where N
        return MultirowMidEnvironment{N}(BC, O, AC)
    end
    function MultirowMidEnvironment(AC::AbstractVector{<:MPSTensor}, O::AbstractMatrix{MPOTensor}, BC::AbstractVector{<:AdjointMPSTensor})
        return MultirowMidEnvironment(BC, O, AC)
    end

    function MultirowMidEnvironment{N}(BC::AbstractVector{<:AdjointMPSTensor}, AC::AbstractVector{<:MPSTensor}, O::AbstractMatrix{MPOTensor}) where N
        return MultirowMidEnvironment{N}(BC, O, AC)
    end
    function MultirowMidEnvironment(BC::AbstractVector{<:AdjointMPSTensor}, AC::AbstractVector{<:MPSTensor}, O::AbstractMatrix{MPOTensor})
        return MultirowMidEnvironment(BC, O, AC)
    end

    function MultirowMidEnvironment{N}(AC::AbstractVector{<:MPSTensor}, BC::AbstractVector{<:AdjointMPSTensor}, O::AbstractMatrix{MPOTensor}) where N
        return MultirowMidEnvironment{N}(BC, O, AC)
    end
    function MultirowMidEnvironment(AC::AbstractVector{<:MPSTensor}, BC::AbstractVector{<:AdjointMPSTensor}, O::AbstractMatrix{MPOTensor})
        return MultirowMidEnvironment(BC, O, AC)
    end

    function MultirowMidEnvironment{N}(O::AbstractMatrix{MPOTensor}, BC::AbstractVector{<:AdjointMPSTensor}, AC::AbstractVector{<:MPSTensor}) where N
        return MultirowMidEnvironment{N}(BC, O, AC)
    end
    function MultirowMidEnvironment(O::AbstractMatrix{MPOTensor}, BC::AbstractVector{<:AdjointMPSTensor}, AC::AbstractVector{<:MPSTensor})
        return MultirowMidEnvironment(BC, O, AC)
    end

    function MultirowMidEnvironment{N}(O::AbstractMatrix{MPOTensor}, AC::AbstractVector{<:MPSTensor}, BC::AbstractVector{<:AdjointMPSTensor}) where N
        return MultirowMidEnvironment{N}(BC, O, AC)
    end
    function MultirowMidEnvironment(O::AbstractMatrix{MPOTensor}, AC::AbstractVector{<:MPSTensor}, BC::AbstractVector{<:AdjointMPSTensor})
        return MultirowMidEnvironment(BC, O, AC)
    end
end
