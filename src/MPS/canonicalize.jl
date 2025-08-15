"""
canonicalize!(obj::DenseInfiniteMPS, si::Int64)

Canonicalize the iMPS where all sites < `si` are left-canonical and all sites > `si` are right-canonical.

`alg::EigenAlgorithm` is the algorithm solving the fixed point.
"""
function canonicalize!(obj::DenseInfiniteMPS{L}, si::Int64; alg::EigenAlgorithm=Defaults.alg_eig) where L
    1 ≤ si ≤ L || throw(ArgumentError("Canonical center position `$si` is out of iMPS range `1 ~ $L`"))

    #! Dev memo: allow default X₀
    _, XL, AL, _ = leftFixedPoint(obj.A, alg)
    _, XR, AR, _ = rightFixedPoint(obj.A, alg)

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
