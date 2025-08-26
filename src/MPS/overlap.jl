"""
    overlap(obj::DenseInfiniteMPS{L, T}, obj′::DenseInfiniteMPS{L, T};
            XL₀::AbstractVector{<:LeftEnvironmentTensor{2}}=_default_X₀_leftFixedPoint(TransferMatrix(getA(obj), adjoint.(getA(obj′)))),
            alg::EigenAlgorithm=Defaults.alg_eig) where {L, T}

Overlap of two iMPSs, i.e., `∏_l λ[l]` solved by
     __                                            __
    |  | -- A[1] -- ... -- A[L] --                |  | --
    |FL|    |              |        =  ∏_l λ[l] * |FL|
    |  | -- B[1] -- ... -- B[L] --                |  | --
     ‾‾                                            ‾‾
Here `A` is `obj.A` or `obj.AL`, and `B` is `adjoint.(obj′.A)` or `adjoint.(obj′.AL)`.
`obj` and `obj′` are normalized iMPS/iMPO.

`abs(overlap(obj, obj′)) = 1` means the two are equivalent in the sense of gauge degrees of freedom.
"""
function overlap(obj::DenseInfiniteMPS{L, T}, obj′::DenseInfiniteMPS{L, T};
            XL₀::AbstractVector{<:LeftEnvironmentTensor{2}}=_default_X₀_leftFixedPoint(TransferMatrix(getAL(obj), adjoint.(getAL(obj′)))),
            alg::EigenAlgorithm=Defaults.alg_eig) where {L, T}
    A = isa(obj, DenseUniformMPS) ? normalize(obj).A : obj.AL
    B = isa(obj′, DenseUniformMPS) ? adjoint.(normalize(obj′).A) : adjoint.(obj′.AL)
    t = TransferMatrix(A, B)
    λ, _, _ = leftFixedPoint(t, XL₀, alg)
    return prod(λ)
end
