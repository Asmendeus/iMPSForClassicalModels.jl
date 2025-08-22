"""
    overlap(obj::T, obj′::T;
            XL₀::AbstractVector{<:LeftEnvironmentTensor{2}}=_default_X₀_leftFixedPoint(TransferMatrix(obj.A, adjoint.(obj′.A))),
            XR₀::AbstractVector{<:RightEnvironmentTensor{2}}=_default_X₀_rightFixedPoint(TransferMatrix(obj.A, adjoint.(obj′.A))),
            alg::EigenAlgorithm=Defaults.alg_eig) where T <: DenseUniformMPS
    overlap(obj::T, obj′::T) where T <: DenseCanonicalMPS

Overlap of two iMPSs or iMPOs, i.e.,
    O(ψ, ψ′) = |⟨ψ|ψ′⟩ / sqrt(⟨ψ|ψ⟩ ⟨ψ′|ψ′⟩)|           (iMPSs)
    O(ρ, ρ′) = |Tr(ρ†ρ′) / sqrt(Tr(ρ†ρ) Tr(ρ′†ρ′))|   (iMPOs)
indicating approximation degree between the two.

`O = 1` means the two are equivalent in the sense of gauge degrees of freedom or one is a submanifold of the other.
"""
function overlap(obj::T, obj′::T;
            XL₀::AbstractVector{<:LeftEnvironmentTensor{2}}=_default_X₀_leftFixedPoint(TransferMatrix(obj.A, adjoint.(obj′.A))),
            XR₀::AbstractVector{<:RightEnvironmentTensor{2}}=_default_X₀_rightFixedPoint(TransferMatrix(obj.A, adjoint.(obj′.A))),
            alg::EigenAlgorithm=Defaults.alg_eig) where T <: DenseUniformMPS
    t = TransferMatrix(obj.A, adjoint.(obj′.A))
    _, XL, _ = leftFixedPoint(t, XL₀, alg)
    _, XR, _ = rightFixedPoint(t, XR₀, alg)
    return norm(contract(environment(XL[1], XR[end]))) / sqrt(norm(obj) * norm(obj′))
end
function overlap(obj::T, obj′::T) where T <: DenseCanonicalMPS
    return tr(obj.C[1].A * obj′.C[1].A') / sqrt(norm(obj) * norm(obj′))
end
