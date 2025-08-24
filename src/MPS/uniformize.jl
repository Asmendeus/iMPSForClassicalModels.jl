"""
    uniformize(obj::DenseCanonicalMPS) -> ::DenseUniformMPS

Uniformize an iMPS with canonical form to an iMPS with uniform form.

# Equations
A[i] = inv(sqrt(C[i-1])) * AL[i] * sqrt(C[i])       # the equation for the application of `uniformize`
     = inv(sqrt(C[i-1])) * AC[i] * inv(sqrt(C[i]))
     = sqrt(C[i-1]) * AR[i] * inv(sqrt(C[i]))
"""
function uniformize(obj::DenseCanonicalMPS{L, T}) where {L, T}
    #> A[i] = inv(sqrt(C[i-1])) * AL[i] * sqrt(C[i])
    right_sqrt_C = sqrt.(obj.C)
    AL = obj.AL
    left_inv_sqrt_C = inv.(right_sqrt_C)[[L, (1:L-1)...]]
    A = left_inv_sqrt_C .* AL .* right_sqrt_C

    #> A[i] = sqrt(C[i-1]) * AR[i] * inv(sqrt(C[i]))
    # right_sqrt_C = sqrt.(obj.C)
    # AR = obj.AR
    # left_sqrt_C = right_sqrt_C[[L, (1:L-1)...]]
    # right_inv_sqrt_C = inv.(right_sqrt_C)
    # A = left_sqrt_C .* AR .* right_inv_sqrt_C
    return uniformize(typeof(obj))(A)
end
uniformize(obj::DenseUniformMPS) = obj
