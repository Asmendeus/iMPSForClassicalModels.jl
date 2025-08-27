"""
    canonicalize(A::AbstractVector{<:AbstractLocalTensor{R}};
                XL₀=_default_X₀_leftFixedPoint(A),
                XR₀=_default_X₀_rightFixedPoint(A),
                alg::EigenAlgorithm=Defaults.alg_eig,
                kwargs...) -> ::InfiniteMPS

Canonicalize a vector of `LocalTensor{R}` (an iMPS/iMPO with uniform form) to an iMPS/iMPO with canonical form

# Steps
1. Solve `AL` and `L` by maximum eigenequation method;
2. Solve `AR` and `R` by maximum eigenequation method;
3. Generate `C` by `C = LR`
4. Generate `AC` by `AL[i] * C[i] = AC[i] = C[i-1] * AR[i]`;
5. Perform a gauge transformation via SVD decomposition to diagonalize `C` and update `AL`, `AR`, `AC`;
"""
function canonicalize(A::AbstractVector{<:AbstractLocalTensor{R}};
                      XL₀=_default_X₀_leftFixedPoint(A),
                      XR₀=_default_X₀_rightFixedPoint(A),
                      alg::EigenAlgorithm=Defaults.alg_eig,
                      kwargs...) where R

    if R == 3
        mpstype = InfiniteMPS
    elseif R == 4
        mpstype = InfiniteMPO
    else
        throw(ArgumentError("Can't convert the vector of `LocalTensor{$R}` to an iMPS(`LocalTensor{3}`) or iMPO(`LocalTensor{4}`)"))
    end
    L = length(A)

    _, XL, AL, _ = leftFixedPoint(A, XL₀, alg; kwargs...)
    _, XR, AR, _ = rightFixedPoint(A, XR₀, alg; kwargs...)

    XL = XL[[(2:L)..., 1]]
    C = XL .* XR

    U = typeof(C)(undef, L)
    Vd = typeof(C)(undef, L)
    for i in 1:L
        u, s, vd = tsvd(C[i].A, (1,), (2,))
        U[i] = eltype(C)(u)
        C[i] = eltype(C)(normalize(s))
        Vd[i] = eltype(C)(vd)
    end
    AL = adjoint.(U)[[L, (1:L-1)...]] .* AL .* U
    AR = Vd[[L, (1:L-1)...]] .* AR .* adjoint.(Vd)

    AC = AL .* C

    return mpstype(AL, AR, AC, C)
end
