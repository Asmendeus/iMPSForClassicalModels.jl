"""
    mutable struct MultirowLeftEnvironment{N} <: AbstractEnvironment{N}
        const BL::AbstractVector{<:AdjointLeftIsometricTensor}
        const O::AbstractMatrix{MPOTensor}
        const AL::AbstractVector{<:LeftIsometricTensor}
    end

Wrapper type for left environment tensor's generating environment.

Graphic presentation:

    —— BL ——
       |
    —— OW ——    (W = N-2)
       |
       ⋮
       |
    —— O1 ——
       |
    —— AL ——

# Constructors
    MultirowLeftEnvironment{N}(BL::AbstractVector{<:AdjointLeftIsometricTensor}, O::AbstractMatrix{MPOTensor}, AL::AbstractVector{<:LeftIsometricTensor})
    MultirowLeftEnvironment(BL::AbstractVector{<:AdjointLeftIsometricTensor}, O::AbstractMatrix{MPOTensor}, AL::AbstractVector{<:LeftIsometricTensor})
    MultirowLeftEnvironment{N}(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractMatrix{MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
    MultirowLeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractMatrix{MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
    MultirowLeftEnvironment{N}(BL::AbstractVector{<:AdjointLeftIsometricTensor}, AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractMatrix{MPOTensor})
    MultirowLeftEnvironment(BL::AbstractVector{<:AdjointLeftIsometricTensor}, AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractMatrix{MPOTensor})
    MultirowLeftEnvironment{N}(AL::AbstractVector{<:LeftIsometricTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}, O::AbstractMatrix{MPOTensor})
    MultirowLeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}, O::AbstractMatrix{MPOTensor})
    MultirowLeftEnvironment{N}(O::AbstractMatrix{MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}, AL::AbstractVector{<:LeftIsometricTensor})
    MultirowLeftEnvironment(O::AbstractMatrix{MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}, AL::AbstractVector{<:LeftIsometricTensor})
    MultirowLeftEnvironment{N}(O::AbstractMatrix{MPOTensor}, AL::AbstractVector{<:LeftIsometricTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
    MultirowLeftEnvironment(O::AbstractMatrix{MPOTensor}, AL::AbstractVector{<:LeftIsometricTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
"""
mutable struct MultirowLeftEnvironment{N} <: AbstractEnvironment{N}
    const BL::AbstractVector{<:AdjointLeftIsometricTensor}
    const O::AbstractMatrix{MPOTensor}
    const AL::AbstractVector{<:LeftIsometricTensor}

    function MultirowLeftEnvironment{N}(BL::AbstractVector{<:AdjointLeftIsometricTensor}, O::AbstractMatrix{MPOTensor}, AL::AbstractVector{<:LeftIsometricTensor}) where N
        length(O) == N - 2 || throw(ArgumentError("The width of multirow MPO doesn't match width of `MultirowLeftEnvironment`"))
        return new{N}(BL, O, AL)
    end
    function MultirowLeftEnvironment(BL::AbstractVector{<:AdjointLeftIsometricTensor}, O::AbstractMatrix{MPOTensor}, AL::AbstractVector{<:LeftIsometricTensor})
        N = length(O) + 2
        return new{N}(BL, O, AL)
    end

    function MultirowLeftEnvironment{N}(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractMatrix{MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}) where N
        return MultirowLeftEnvironment{N}(BL, O, AL)
    end
    function MultirowLeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractMatrix{MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
        return MultirowLeftEnvironment(BL, O, AL)
    end

    function MultirowLeftEnvironment{N}(BL::AbstractVector{<:AdjointLeftIsometricTensor}, AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractMatrix{MPOTensor}) where N
        return MultirowLeftEnvironment{N}(BL, O, AL)
    end
    function MultirowLeftEnvironment(BL::AbstractVector{<:AdjointLeftIsometricTensor}, AL::AbstractVector{<:LeftIsometricTensor}, O::AbstractMatrix{MPOTensor})
        return MultirowLeftEnvironment(BL, O, AL)
    end

    function MultirowLeftEnvironment{N}(AL::AbstractVector{<:LeftIsometricTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}, O::AbstractMatrix{MPOTensor}) where N
        return MultirowLeftEnvironment{N}(BL, O, AL)
    end
    function MultirowLeftEnvironment(AL::AbstractVector{<:LeftIsometricTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}, O::AbstractMatrix{MPOTensor})
        return MultirowLeftEnvironment(BL, O, AL)
    end

    function MultirowLeftEnvironment{N}(O::AbstractMatrix{MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}, AL::AbstractVector{<:LeftIsometricTensor}) where N
        return MultirowLeftEnvironment{N}(BL, O, AL)
    end
    function MultirowLeftEnvironment(O::AbstractMatrix{MPOTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}, AL::AbstractVector{<:LeftIsometricTensor})
        return MultirowLeftEnvironment(BL, O, AL)
    end

    function MultirowLeftEnvironment{N}(O::AbstractMatrix{MPOTensor}, AL::AbstractVector{<:LeftIsometricTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor}) where N
        return MultirowLeftEnvironment{N}(BL, O, AL)
    end
    function MultirowLeftEnvironment(O::AbstractMatrix{MPOTensor}, AL::AbstractVector{<:LeftIsometricTensor}, BL::AbstractVector{<:AdjointLeftIsometricTensor})
        return MultirowLeftEnvironment(BL, O, AL)
    end
end
