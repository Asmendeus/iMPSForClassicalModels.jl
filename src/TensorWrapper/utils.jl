"""
    function isAdjoint(::AbstractTensorWrapperper)::Bool

Whether the tensor wrapper is an adjoint wrapper
"""
isAdjoint(::Union{LocalTensor, LeftIsometricTensor, RightIsometricTensor, LeftEnvironmentTensor, RightEnvironmentTensor}) = false
isAdjoint(::Union{AdjointLocalTensor, AdjointLeftIsometricTensor, AdjointRightIsometricTensor}) = true

"""
    function isLeftIsometric(A::AbstractTensorMap, isadjoint::Bool; tol::Float64=1e-8)
    function isLeftIsometric(A::AbstractLocalTensor; tol::Float64=1e-8)

# Return
-`::Bool`: whether the R-leg local tensor `A` is left orthogonal
"""
function isLeftIsometric(A::AbstractTensorMap, isadjoint::Bool; tol::Float64=Defaults.tol_norm)
    if numout(A) == 1 && numin(A) == 1 && isadjoint == false
        # Union{BondTensor, LeftIsometricBondTensor, RightIsometricBondTensor}
        @tensor E[-1; -2] := A[1 -2] * A'[-1 1]
        return norm(E - id(domain(A))) < tol
    elseif numout(A) == 1 && numin(A) == 1 && isadjoint == true
        # Union{AdjointBondTensor, AdjointLeftIsometricBondTensor, AdjointRightIsometricBondTensor}
        @tensor E[-1; -2] := A'[1 -2] * A[-1 1]
        return norm(E - id(codomain(A))) < tol
    elseif numout(A) == 2 && numin(A) == 1 && isadjoint == false
        # Union{MPSTensor, LeftIsometricMPSTensor, RightIsometricMPSTensor}
        @tensor E[-1; -2] := A[1 2 -2] * A'[-1 1 2]
        return norm(E - id(domain(A))) < tol
    elseif numout(A) == 1 && numin(A) == 2 && isadjoint == true
        # Union{AdjointMPSTensor, AdjointLeftIsometricMPSTensor, AdjointRightIsometricMPSTensor}
        @tensor E[-1; -2] := A'[1 2 -2] * A[-1 1 2]
        return norm(E - id(codomain(A))) < tol
    elseif numout(A) == 2 && numin(A) == 2 && isadjoint == false
        # Union{MPOTensor, LeftIsometricMPOTensor, RightIsometricMPOTensor}
        @tensor E[-1; -2] := A[1 2 3 -2] * A'[3 -1 1 2]
        return norm(E - id(domain(A, 2))) < tol
    elseif numout(A) == 2 && numin(A) == 2 && isadjoint == true
        # Union{AdjointMPOTensor, AdjointLeftIsometricMPOTensor, AdjointRightIsometricMPOTensor}
        @tensor E[-1; -2] := A'[1 2 3 -2] * A[3 -1 1 2]
        return norm(E - id(codomain(A, 2))) < tol
    else
        # R > 4
        throw(ArgumentError("Behavior not yet defined"))
    end
end
function isLeftIsometric(A::AbstractLocalTensor; tol::Float64=Defaults.tol_norm)::Bool
    return isLeftIsometric(A.A, isAdjoint(A); tol=tol)
end


"""
    function isRightIsometric(A::AbstractTensorMap, isadjoint::Bool; tol::Float64=1e-8)
    function isRightIsometric(A::AbstractLocalTensor; tol::Float64=1e-8)

# Return
-`::Bool`: whether the R-leg local tensor `A` is right orthogonal
"""
function isRightIsometric(A::AbstractTensorMap, isadjoint::Bool; tol::Float64=1e-8)
    if numout(A) == 1 && numin(A) == 1 && isadjoint == false
        # Union{BondTensor, LeftIsometricBondTensor, RightIsometricBondTensor}
        @tensor E[-1; -2] := A[-1 1] * A'[1 -2]
        return norm(E - id(codomain(A, 1))) < tol
    elseif numout(A) == 1 && numin(A) == 1 && isadjoint == true
        # Union{AdjointBondTensor, AdjointLeftIsometricBondTensor, AdjointRightIsometricBondTensor}
        @tensor E[-1; -2] := A'[-1 1] * A[1 -2]
        return norm(E - id(domain(A, 1))) < tol
    elseif numout(A) == 2 && numin(A) == 1 && isadjoint == false
        # Union{MPSTensor, LeftIsometricMPSTensor, RightIsometricMPSTensor}
        @tensor E[-1; -2] := A[-1 1 2] * A'[2 -2 1]
        return norm(E - id(codomain(A, 1))) < tol
    elseif numout(A) == 1 && numin(A) == 2 && isadjoint == true
        # Union{AdjointMPSTensor, AdjointLeftIsometricMPSTensor, AdjointRightIsometricMPSTensor}
        @tensor E[-1; -2] := A'[-1 1 2] * A[2 -2 1]
        return norm(E - id(domain(A, 1))) < tol
    elseif numout(A) == 2 && numin(A) == 2 && isadjoint == false
        # Union{MPOTensor, LeftIsometricMPOTensor, RightIsometricMPOTensor}
        @tensor E[-1; -2] := A[-1 1 2 3] * A'[2 3 -2 1]
        return norm(E - id(codomain(A, 1))) < tol
    elseif numout(A) == 2 && numin(A) == 2 && isadjoint == true
        # Union{AdjointMPOTensor, AdjointLeftIsometricMPOTensor, AdjointRightIsometricMPOTensor}
        @tensor E[-1; -2] := A'[-1 1 2 3] * A[2 3 -2 1]
        return norm(E - id(domain(A, 1))) < tol
    else
        # R > 4
        throw(ArgumentError("Behavior not yet defined"))
    end
end
function isRightIsometric(A::AbstractLocalTensor; tol::Float64=Defaults.tol_norm)::Bool
    return isRightIsometric(A.A, isAdjoint(A); tol=tol)
end

"""
     leftorth(A::LocalTensor{R₁}; trunc = notrunc(), kwargs...)
     leftorth(A::AdjointLocalTensor{R₁}; trunc = notrunc(), kwargs...)

Left canonicalize a local tensor.

# Return
-`::LeftIsometricTensor{R₁}`: left isometric tensor
-`::BondTensor`: bond tensor

If `trunc = notrunc()`, use `TensorKit.leftorth`, otherwise, use `TensorKit.tsvd`.
Propagate `kwargs` to the TensorKit functions.
"""
function leftorth(A::LocalTensor{R₁}; trunc=notrunc(), kwargs...) where R₁
    if trunc == notrunc()
        Q, R = leftorth(A.A, (Tuple(1:R₁-1), (R₁,)); kwargs...)
        if R₁ == 2 || R₁ == 3
            return LeftIsometricTensor(Q; check=false), LocalTensor(R), BondInfo(Q, :R)
        else
            return LeftIsometricTensor(permute(Q, (1, 2), Tuple(3:R₁)); check=false), LocalTensor(R), BondInfo(Q, :R)
        end
    else
        u, s, vd, info = tsvd(A, (Tuple(1:R₁-1), (R₁,)); trunc=trunc, kwargs...)
        if R₁ == 2 || R₁ == 3
            return LeftIsometricTensor(u; check=false), LocalTensor(s * vd), info
        else
            return LeftIsometricTensor(permute(u, (1, 2), Tuple(3:R₁)); check=false), LocalTensor(s * vd), info
        end
    end
end
function leftorth(A::AdjointLocalTensor{R₁}; trunc=notrunc(), kwargs...) where R₁
    if trunc == notrunc()
        if R₁ == 2
            Q, R = leftorth(A.A, ((2,), (1,)); kwargs...)
            return AdjointLeftIsometricTensor(permute(Q, (2,), (1,)); check=false), AdjointLocalTensor(permute(R, (2,), (1,))), BondInfo(Q, :R)
        else
            Q, R = leftorth(A.A, (Tuple(vcat(1:R₁-3, [R₁-1, R₁])), (R₁-2,)); kwargs...)
            return AdjointLeftIsometricTensor(permute(Q, Tuple(vcat(1:R₁-3, [R₁,])), (R₁-2, R₁-1)); check=false), AdjointLocalTensor(permute(R, (2,), (1,))), BondInfo(Q, :R)
        end
    else
        if R₁ == 2
            u, s, vd, info = tsvd(A, ((2,), (1,)); trunc=trunc, kwargs...)
            return AdjointLeftIsometricTensor(permute(u, (2,), (1,)); check=false), AdjointLocalTensor(permute(s * vd, (2,), (1,))), info
        else
            u, s, vd, info =  tsvd(A, (Tuple(vcat(1:R₁-3, [R₁-1, R₁])), (R₁-2,)); trunc=trunc, kwargs...)
            return AdjointLeftIsometricTensor(permute(u, Tuple(vcat(1:R₁-3, [R₁,])), (R₁-2, R₁-1)); check=false), AdjointLocalTensor(permute(s * vd, (2,), (1,))), info
        end
    end
end

"""
    rightorth{R₂}(A::LocalTensor; trunc = notrunc(), kwargs...)
    rightorth{R₂}(A::AdjointLocalTensor; trunc = notrunc(), kwargs...)

Right canonicalize a local tensor.

# Return
-`::BondTensor`: bond tensor
-`::RightIsometricTensor{R₂}`: right isometric tensor

If `trunc = notrunc()`, use `TensorKit.rightorth`, otherwise, use `TensorKit.tsvd`.
Propagate `kwargs` to the TensorKit functions.
"""
function rightorth(A::LocalTensor{R₂}; trunc=notrunc(), kwargs...) where R₂
    if trunc == notrunc()
        L, Q = rightorth(A.A, ((1,), Tuple(2:R₂)); kwargs...)
        if R₂ == 2
            return LocalTensor(L), RightIsometricTensor(Q; check=false), BondInfo(Q, :L)
        else
            return LocalTensor(L), RightIsometricTensor(permute(Q, (1, 2), Tuple(3:R₂)); check=false), BondInfo(Q, :L)
        end
    else
        u, s, vd, info = tsvd(A, ((1,), Tuple(2:R₂)); trunc=trunc, kwargs...)
        if R₂ == 2
            return LocalTensor(u * s), RightIsometricTensor(vd; check=false), info
        else
            return LocalTensor(u * s), RightIsometricTensor(permute(vd, (1, 2), Tuple(3:R₂)); check=false), info
        end
    end
end
function rightorth(A::AdjointLocalTensor{R₂}; trunc=notrunc(), kwargs...) where R₂
    if trunc == notrunc()
        if R₂ == 2
            L, Q = rightorth(A.A, ((2,), (1, )); kwargs...)
            return LocalTensor(permute(L, (2,), (1,))), AdjointRightIsometricTensor(permute(Q, (2,), (1,)); check=false), BondInfo(Q, :L)
        else
            L, Q = rightorth(A.A, ((R₂-1,), Tuple(vcat((1:R₂-2), [R₂,]))); kwargs...)
            return LocalTensor(permute(L, (2,), (1,))), AdjointRightIsometricTensor(permute(Q, Tuple(2:R₂-1), (1, R₂)); check=false), BondInfo(Q, :L)
        end
    else
        if R₂ == 2
            u, s, vd, info = tsvd(A, ((2,), (1, )); trunc=trunc, kwargs...)
            return LocalTensor(permute(u * s, (2,), (1,))), AdjointRightIsometricTensor(permute(vd, (2,), (1,)); check=false), info
        else
            u, s, vd, info = tsvd(A, ((R₂-1,), Tuple(vcat((1:R₂-2), [R₂,]))); trunc=trunc, kwargs...)
            return LocalTensor(permute(u * s, (2,), (1,))), AdjointRightIsometricTensor(permute(vd, Tuple(2:R₂-1), (1, R₂)); check=false), info
        end
    end
end
