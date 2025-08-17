"""
    canonicalize!(obj::DenseInfiniteMPS{L, T}, si::Int64; XL₀::AbstractVector{<:LocalTensor{R}}=_default_X₀_leftFixedPoint(obj.A), XR₀::AbstractVector{<:LocalTensor{R}}=_default_X₀_rightFixedPoint(obj.A), alg::EigenAlgorithm=Defaults.alg_eig)

Canonicalize an iMPS/iMPO with uniform form.

`alg::EigenAlgorithm` is the algorithm solving the fixed point.
"""
function canonicalize!(obj::DenseInfiniteMPS{L, T}, si::Int64; XL₀::AbstractVector{<:LocalTensor{R}}=_default_X₀_leftFixedPoint(obj.A), XR₀::AbstractVector{<:LocalTensor{R}}=_default_X₀_rightFixedPoint(obj.A), alg::EigenAlgorithm=Defaults.alg_eig) where {L, T, R}

    1 ≤ si ≤ L || throw(ArgumentError("Canonical center position `$si` is out of iMPS range `1 ~ $L`"))
    isuniform(obj) || throw(ArgumentError("Illegal behavior of canonicalizing an iMPS/iMPO with canonical form"))

    _, XL, AL, _ = leftFixedPoint(obj.A, XL₀, alg)
    _, XR, AR, _ = rightFixedPoint(obj.A, XR₀, alg)

    for i in 1:si-1
        obj.A[i] = AL[i]
    end

    for i in si+1:L
        obj.A[i] = AR[i]
    end

    obj.A[si] = XL[si] * obj.A[si] * XR[si]
    normalize!(obj.A[si])

    obj.Center = si
    return obj
end

function canonicalize(obj::DenseInfiniteMPS{L, T}, si::Int64; XL₀::AbstractVector{<:LocalTensor{R}}=_default_X₀_leftFixedPoint(obj.A), XR₀::AbstractVector{<:LocalTensor{R}}=_default_X₀_rightFixedPoint(obj.A), alg::EigenAlgorithm=Defaults.alg_eig) where {L, T, R}
    obj′ = deepcopy(obj)
    return canonicalize!(obj′, si; XL₀=XL₀, XR₀=XR₀, alg=alg)
end
