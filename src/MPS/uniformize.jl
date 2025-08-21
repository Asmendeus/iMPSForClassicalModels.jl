"""
    uniformize(obj::DenseCanonicalMPS) -> ::DenseUniformMPS

Uniformize an iMPS with canonical form to an iMPS with uniform form
"""
function uniformize(obj::DenseCanonicalMPS)
    return uniformize(typeof(obj))(A)
end
