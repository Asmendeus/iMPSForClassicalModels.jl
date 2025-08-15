"""
canonicalize!(obj::DenseInfiniteMPS, si::Int64; kwargs...)

Canonicalize the iMPS where all sites < `si` are left-canonical and all sites > `si` are right-canonical.

`alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}` is the algorithm solving the fixed point.
`kwargs` is propagated to `leftorth` and `rightorth` to determine how to truncate the SVD spectra.
"""
function canonicalize!(obj::DenseInfiniteMPS{L}, si::Int64; alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Defaults.alg_eig(), kwargs...) where L
    1 ≤ si ≤ L || throw(ArgumentError("Canonical center position `$si` is out of iMPS range `1 ~ $L`"))

    _, XL, AL, _ = leftFixedPoint(obj.A, alg)
    _, XR, AR, _ = rightFixedPoint(obj.A, alg)

    for i in 1:si-1
        obj.A[i] = AR[i]
    end

    for i in si+1:L
        obj.A[i] = AL[i]
    end

    obj.A[si] = XL[si] * obj.A[si] * XR[si]
    normalize!(obj.A[si])

    Center(obj)[:] = [si, si]
    return obj
end
