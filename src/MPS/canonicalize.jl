"""
    canonicalize(obj::DenseUniformMPS{L, T}; XL₀=_default_X₀_leftFixedPoint(obj.A),
            XR₀=_default_X₀_rightFixedPoint(obj.A), alg::EigenAlgorithm=Defaults.alg_eig,
            kwargs...) -> ::DenseCanonicalMPS

Canonicalize an iMPS with uniform form to an iMPS with canonical form

# Steps
1. Solve `AL` and `L` by maximum eigenequation method;
2. Solve `AR` and `R` by maximum eigenequation method;
3. Generate `C` by `C = LR`
4. Generate `AC` by `AL[i] * C[i] = AC[i] = C[i-1] * AR[i]`;
5. Perform a gauge transformation via SVD decomposition to diagonalize `C` and update `AL`, `AR`, `AC`;
"""
function canonicalize(obj::DenseUniformMPS{L, T}; XL₀=_default_X₀_leftFixedPoint(obj.A),
            XR₀=_default_X₀_rightFixedPoint(obj.A), alg::EigenAlgorithm=Defaults.alg_eig, kwargs...) where {L, T}

    @assert alg isa KrylovKit.KrylovAlgorithm "Only the maximum eigenequation method is supported to canonicalize."

    _, XL, AL, _ = leftFixedPoint(obj.A, XL₀, alg; kwargs...)
    _, XR, AR, _ = rightFixedPoint(obj.A, XR₀, alg; kwargs...)

    XL = XL[[(2:L)..., 1]]
    C = XL .* XR

    U = typeof(C)(undef, L)
    Vd = typeof(C)(undef, L)
    for i in 1:L
        u, s, vd = tsvd(C[i].A, (1,), (2,))
        U[i] = eltype(C)(u)
        C[i] = eltype(C)(s)
        Vd[i] = eltype(C)(vd)
    end
    AL = adjoint.(U)[[L, (1:L-1)...]] .* AL .* U
    AR = Vd[[L, (1:L-1)...]] .* AR .* adjoint.(Vd)

    AC = AL .* C

    return canonicalize(typeof(obj))(AL, AR, AC, C)
end
canonicalize(obj::DenseCanonicalMPS; kwargs...) = obj
