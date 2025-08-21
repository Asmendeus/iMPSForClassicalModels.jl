"""
    canonicalize(obj::DenseUniformMPS) -> ::DenseCanonicalMPS

Canonicalize an iMPS with uniform form to an iMPS with canonical form
"""
function canonicalize(obj::DenseUniformMPS)
    return canonicalize(typeof(obj))(AL, AR, AC, C)
end
