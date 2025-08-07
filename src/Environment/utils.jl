function Environment(B::AbstractVector{<:AbstractTensorWrapper}, O::AbstractVector{<:AbstractTensorWrapper}, A::AbstractVector{<:AbstractTensorWrapper})
    if Union{eltype(A), eltype(B), eltype(O)} <: Union{LeftIsometricMPSTensor, AdjointLeftIsometricMPSTensor, MPOTensor}
        return LeftEnvironment(B, O, A)
    elseif Union{eltype(A), eltype(B), eltype(O)} <: Union{RightIsometricMPSTensor, AdjointRightIsometricMPSTensor, MPOTensor}
        return RightEnvironment(B, O, A)
    else
        throw(ArgumentError("unsupported combinations of types (::$(typeof(B)), ::$(typeof(O))), ::$(typeof(A))"))
    end
end

function Environment(B::AbstractVector{<:AbstractTensorWrapper}, O::AbstractMatrix{<:AbstractTensorWrapper}, A::AbstractVector{<:AbstractTensorWrapper})
    if Union{eltype(A), eltype(B), eltype(O)} <: Union{LeftIsometricMPSTensor, AdjointLeftIsometricMPSTensor, MPOTensor}
        return MultiLeftEnvironment(B, O, A)
    elseif Union{eltype(A), eltype(B), eltype(O)} <: Union{RightIsometricMPSTensor, AdjointRightIsometricMPSTensor, MPOTensor}
        return MultiRightEnvironment(B, O, A)
    else
        throw(ArgumentError("unsupported combinations of types (::$(typeof(B)), ::$(typeof(O))), ::$(typeof(A))"))
    end
end
