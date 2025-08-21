"""
    canonicalize(obj::DenseUniformMPS{L, T}; XL₀=_default_X₀_leftFixedPoint(obj.A),
            XR₀=_default_X₀_rightFixedPoint(obj.A), alg::EigenAlgorithm=Defaults.alg_eig,
            kwargs...) -> ::DenseCanonicalMPS

Canonicalize an iMPS with uniform form to an iMPS with canonical form

# Steps
1. Solve `AL` and `L` by maximum eigenequation method;
2. Solve `AR` and `R` by maximum eigenequation method;
3. Generate `C` by `C = LR`
4. Perform a gauge transformation via SVD decomposition to diagonalize `C` and update `AL` and `AR`;
5. Generate `AC` by `AL[i] * C[i] = AC[i] = C[i-1] * AR[i]`.
"""
function canonicalize(obj::DenseUniformMPS{L, T}; XL₀=_default_X₀_leftFixedPoint(obj.A),
            XR₀=_default_X₀_rightFixedPoint(obj.A), alg::EigenAlgorithm=Defaults.alg_eig, kwargs...) where {L, T}

    @assert alg isa KrylovKit.KrylovAlgorithm "Only the maximum eigenequation method is supported to canonicalize."

    _, XL, AL, _ = leftFixedPoint(obj.A, XL₀, alg; kwargs...)
    _, XR, AR, _ = rightFixedPoint(obj.A, XR₀, alg; kwargs...)

    XL = XL[(2:end)..., 1]
    C = XL .* XR

    U = Vector{typeof(C)}(undef, L)
    Vd = Vector{typeof(C)}(undef, L)
    for i in 1:L
        u, s, vd = tsvd(C[i].A, (1,), (2,))
        U[i] = typeof(C)(u)
        C[i] = typeof(C)(s)
        Vd[i] = typeof(C)(vd)
    end
    AL = AL .* U
    AR = Vd[end, (1:end-1)...] .* AR

    AC = AL .* C

    return canonicalize(typeof(obj))(AL, AR, AC, C)
end
