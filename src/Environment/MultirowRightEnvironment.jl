"""
    mutable struct MultirowRightEnvironment{N} <: AbstractEnvironment{N}
        const BR::AbstractVector{<:AdjointRightIsometricTensor}
        const O::AbstractMatrix{MPOTensor}
        const AR::AbstractVector{<:RightIsometricTensor}
    end

Wrapper type for Right environment tensor's generating environment.

Graphic presentation:

    —— BR ——
       |
    —— OW ——    (W = N-2)
       |
       ⋮
       |
    —— O1 ——
       |
    —— AR ——

# Constructors
    MultirowRightEnvironment{N}(BR::AbstractVector{<:AdjointRightIsometricTensor}, O::AbstractMatrix{MPOTensor}, AR::AbstractVector{<:RightIsometricTensor})
    MultirowRightEnvironment(BR::AbstractVector{<:AdjointRightIsometricTensor}, O::AbstractMatrix{MPOTensor}, AR::AbstractVector{<:RightIsometricTensor})
    MultirowRightEnvironment{N}(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractMatrix{MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
    MultirowRightEnvironment(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractMatrix{MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
    MultirowRightEnvironment{N}(BR::AbstractVector{<:AdjointRightIsometricTensor}, AR::AbstractVector{<:RightIsometricTensor}, O::AbstractMatrix{MPOTensor})
    MultirowRightEnvironment(BR::AbstractVector{<:AdjointRightIsometricTensor}, AR::AbstractVector{<:RightIsometricTensor}, O::AbstractMatrix{MPOTensor})
    MultirowRightEnvironment{N}(AR::AbstractVector{<:RightIsometricTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}, O::AbstractMatrix{MPOTensor})
    MultirowRightEnvironment(AR::AbstractVector{<:RightIsometricTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}, O::AbstractMatrix{MPOTensor})
    MultirowRightEnvironment{N}(O::AbstractMatrix{MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}, AR::AbstractVector{<:RightIsometricTensor})
    MultirowRightEnvironment(O::AbstractMatrix{MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}, AR::AbstractVector{<:RightIsometricTensor})
    MultirowRightEnvironment{N}(O::AbstractMatrix{MPOTensor}, AR::AbstractVector{<:RightIsometricTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
    MultirowRightEnvironment(O::AbstractMatrix{MPOTensor}, AR::AbstractVector{<:RightIsometricTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
"""
mutable struct MultirowRightEnvironment{N} <: AbstractEnvironment{N}
    const BR::AbstractVector{<:AdjointRightIsometricTensor}
    const O::AbstractMatrix{MPOTensor}
    const AR::AbstractVector{<:RightIsometricTensor}

    function MultirowRightEnvironment{N}(BR::AbstractVector{<:AdjointRightIsometricTensor}, O::AbstractMatrix{MPOTensor}, AR::AbstractVector{<:RightIsometricTensor}) where N
        length(O) == N - 2 || throw(ArgumentError("The width of multirow MPO doesn't match width of `MultirowRightEnvironment`"))
        return new{N}(BR, O, AR)
    end
    function MultirowRightEnvironment(BR::AbstractVector{<:AdjointRightIsometricTensor}, O::AbstractMatrix{MPOTensor}, AR::AbstractVector{<:RightIsometricTensor})
        N = length(O) + 2
        return new{N}(BR, O, AR)
    end

    function MultirowRightEnvironment{N}(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractMatrix{MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}) where N
        return MultirowRightEnvironment{N}(BR, O, AR)
    end
    function MultirowRightEnvironment(AR::AbstractVector{<:RightIsometricTensor}, O::AbstractMatrix{MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
        return MultirowRightEnvironment(BR, O, AR)
    end

    function MultirowRightEnvironment{N}(BR::AbstractVector{<:AdjointRightIsometricTensor}, AR::AbstractVector{<:RightIsometricTensor}, O::AbstractMatrix{MPOTensor}) where N
        return MultirowRightEnvironment{N}(BR, O, AR)
    end
    function MultirowRightEnvironment(BR::AbstractVector{<:AdjointRightIsometricTensor}, AR::AbstractVector{<:RightIsometricTensor}, O::AbstractMatrix{MPOTensor})
        return MultirowRightEnvironment(BR, O, AR)
    end

    function MultirowRightEnvironment{N}(AR::AbstractVector{<:RightIsometricTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}, O::AbstractMatrix{MPOTensor}) where N
        return MultirowRightEnvironment{N}(BR, O, AR)
    end
    function MultirowRightEnvironment(AR::AbstractVector{<:RightIsometricTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}, O::AbstractMatrix{MPOTensor})
        return MultirowRightEnvironment(BR, O, AR)
    end

    function MultirowRightEnvironment{N}(O::AbstractMatrix{MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}, AR::AbstractVector{<:RightIsometricTensor}) where N
        return MultirowRightEnvironment{N}(BR, O, AR)
    end
    function MultirowRightEnvironment(O::AbstractMatrix{MPOTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}, AR::AbstractVector{<:RightIsometricTensor})
        return MultirowRightEnvironment(BR, O, AR)
    end

    function MultirowRightEnvironment{N}(O::AbstractMatrix{MPOTensor}, AR::AbstractVector{<:RightIsometricTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor}) where N
        return MultirowRightEnvironment{N}(BR, O, AR)
    end
    function MultirowRightEnvironment(O::AbstractMatrix{MPOTensor}, AR::AbstractVector{<:RightIsometricTensor}, BR::AbstractVector{<:AdjointRightIsometricTensor})
        return MultirowRightEnvironment(BR, O, AR)
    end
end
