"""
    function isLeftIsometric(A::AbstractTensorMap; tol::Float64=1e-8)
    function isLeftIsometric(A::AbstractMPSTensor; tol::Float64=1e-8)

# Return
-`::Bool`: whether the 3-leg MPS tensor `A` is left orthogonal
"""
function isLeftIsometric(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm)::Bool
    if length(codomain(A)) == 2 && length(domain(A)) == 1
        @tensor E[-1; -2] := A[1 2 -1] * A'[-2 1 2]
        return norm(E - id(domain(A))) < tol
    elseif length(codomain(A)) == 1 && length(domain(A)) == 2
        @tensor E[-1; -2] := A'[1 2 -1] * A[-2 1 2]
        return norm(E - id(codomain(A))) < tol
    else
        throw(ArgumentError("The space $(space(A)) does not conform to `AbstractMPSTensor` space"))
    end
end
function isLeftIsometric(A::AbstractMPSTensor; tol::Float64=Defaults.tol_norm)::Bool
    return isLeftIsometric(A.A; tol=tol)
end

"""
    function isRightIsometric(A::AbstractTensorMap; tol::Float64=1e-8)
    function isRightIsometric(A::AbstractMPSTensor; tol::Float64=1e-8)

# Return
-`::Bool`: whether the 3-leg MPS tensor `A` is right orthogonal
"""
function isRightIsometric(A::AbstractTensorMap; tol::Float64=Defaults.tol_norm)::Bool
    if length(codomain(A)) == 2 && length(domain(A)) == 1
        @tensor E[-1; -2] := A[-1 1 2] * A'[2 -2 1]
        return norm(E - id(codomain(A)[1])) < tol
    elseif length(codomain(A)) == 1 && length(domain(A)) == 2
        @tensor E[-1; -2] := A'[-1 1 2] * A[2 -2 1]
        return norm(E - id(domain(A)[1])) < tol
    else
        throw(ArgumentError("The space $(space(A)) does not conform to `AbstractMPSTensor` space"))
    end
end
function isRightIsometric(A::AbstractMPSTensor; tol::Float64=Defaults.tol_norm)::Bool
    return isRightIsometric(A.A; tol=tol)
end

"""
     leftorth(A::MPSTensor; trunc = notrunc(), kwargs...)
     leftorth(A::AdjointMPSTensor; trunc = notrunc(), kwargs...)

Left canonicalize a on-site MPS tensor.

# Return
-`Q::AbstractTensorMap`: left isometric tensor
-`R::AbstractTensorMap`: upper triangular tensor
-`info::BondInfo`: bond information

If `trunc = notrunc()`, use `TensorKit.leftorth`, otherwise, use `TensorKit.tsvd`.
Propagate `kwargs` to the TensorKit functions.
"""
function leftorth(A::MPSTensor; trunc=notrunc(), kwargs...)
     if trunc == notrunc()
          Q, R = leftorth(A.A, ((1, 2), (3,)); kwargs...)
          return LeftIsometricMPSTensor(Q; check=false), R
     else
          u, s, vd, _ =  tsvd(A, ((1, 2), (3,)); trunc=trunc, kwargs...)
          return LeftIsometricMPSTensor(u; check=false), s * vd
     end
end
function leftorth(A::AdjointMPSTensor; trunc=notrunc(), kwargs...)
     if trunc == notrunc()
          Q, R = leftorth(A.A, ((2, 3), (1,)); kwargs...)
          return AdjointLeftIsometricMPSTensor(permute(Q, (3,), (1, 2)); check=false), R
     else
          u, s, vd, _ =  tsvd(A, ((2, 3), (1,)); trunc=trunc, kwargs...)
          return AdjointLeftIsometricMPSTensor(permute(u, (3,), (1, 2)); check=false), s * vd
     end
end

"""
    rightorth(A::MPSTensor; trunc = notrunc(), kwargs...)
    rightorth(A::AdjointMPSTensor; trunc = notrunc(), kwargs...)

Right canonicalize a on-site MPS tensor.

# Return
-`L::AbstractTensorMap`: lower triangular tensor
-`Q::AbstractTensorMap`: right isometric tensor
-`info::BondInfo`: bond information

If `trunc = notrunc()`, use `TensorKit.rightorth`, otherwise, use `TensorKit.tsvd`.
Propagate `kwargs` to the TensorKit functions.
"""
function rightorth(A::MPSTensor; trunc=notrunc(), kwargs...)
     if trunc == notrunc()
          L, Q = rightorth(A.A, ((1,), (2, 3)); kwargs...)
          return L, RightIsometricMPSTensor(permute(Q, (1, 2), (3, )); check=false)
     else
          u, s, vd, _ = tsvd(A, ((1,), (2, 3)); trunc=trunc, kwargs...)
          return u * s, RightIsometricMPSTensor(permute(vd, (1, 2), (3, )); check=false)
     end
end
function rightorth(A::AdjointMPSTensor; trunc=notrunc(), kwargs...)
     if trunc == notrunc()
          L, Q = rightorth(A.A, ((2,), (3, 1)); kwargs...)
          return L, AdjointRightIsometricMPSTensor(permute(Q, (3,), (1, 2)); check=false)
     else
          u, s, vd, _ = tsvd(A, ((2,), (3, 1)); trunc=trunc, kwargs...)
          return u * s, AdjointRightIsometricMPSTensor(permute(vd, (3,), (1, 2)); check=false)
     end
end
