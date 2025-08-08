"""
    function isLeftIsometric(A::AbstractMPSTensor; tol::Float64=1e-8)

# Return
-`::Bool`: whether the R-leg MPS tensor `A` is left orthogonal
"""
function isLeftIsometric(A::AbstractMPSTensor; tol::Float64=Defaults.tol_norm)::Bool
    if A isa Union{MPSTensor{3}, LeftIsometricTensor{3}, RightIsometricTensor{3}}
        @tensor E[-1; -2] := A.A[1 2 -2] * A.A'[-1 1 2]
        return norm(E - id(domain(A))) < tol
    elseif A isa Union{AdjointMPSTensor{3}, AdjointLeftIsometricTensor{3}, AdjointRightIsometricTensor{3}}
        @tensor E[-1; -2] := A.A'[1 2 -2] * A.A[-1 1 2]
        return norm(E - id(codomain(A))) < tol
    elseif A isa Union{MPSTensor{4}, LeftIsometricTensor{4}, RightIsometricTensor{4}}
        @tensor E[-1; -2] := A.A[1 2 3 -2] * A.A'[3 -1 1 2]
        return norm(E - id(domain(A, 2))) < tol
    elseif A isa Union{AdjointMPSTensor{4}, AdjointLeftIsometricTensor{4}, AdjointRightIsometricTensor{4}}
        @tensor E[-1; -2] := A.A'[1 2 3 -2] * A.A[3 -1 1 2]
        return norm(E - id(codomain(A, 2))) < tol
    else
        # R == 2 || R > 4
        throw(ArgumentError("Behavior not yet defined"))
    end
end

"""
    function isRightIsometric(A::AbstractMPSTensor; tol::Float64=1e-8)

# Return
-`::Bool`: whether the R-leg MPS tensor `A` is right orthogonal
"""
function isRightIsometric(A::AbstractMPSTensor; tol::Float64=Defaults.tol_norm)::Bool
    if A isa Union{MPSTensor{3}, LeftIsometricTensor{3}, RightIsometricTensor{3}}
        @tensor E[-1; -2] := A.A[-1 1 2] * A.A'[2 -2 1]
        return norm(E - id(codomain(A, 1))) < tol
    elseif A isa Union{AdjointMPSTensor{3}, AdjointLeftIsometricTensor{3}, AdjointRightIsometricTensor{3}}
        @tensor E[-1; -2] := A.A'[-1 1 2] * A.A[2 -2 1]
        return norm(E - id(domain(A, 1))) < tol
    elseif A isa Union{MPSTensor{4}, LeftIsometricTensor{4}, RightIsometricTensor{4}}
        @tensor E[-1; -2] := A.A[-1 1 2 3] * A.A'[2 3 -2 1]
        return norm(E - id(codomain(A, 1))) < tol
    elseif A isa Union{AdjointMPSTensor{4}, AdjointLeftIsometricTensor{4}, AdjointRightIsometricTensor{4}}
        @tensor E[-1; -2] := A.A'[-1 1 2 3] * A.A[2 3 -2 1]
        return norm(E - id(domain(A, 1))) < tol
    else
        # R == 2 || R > 4
        throw(ArgumentError("Behavior not yet defined"))
    end
end

"""
     leftorth(A::MPSTensor{R₁}; trunc = notrunc(), kwargs...)
     leftorth(A::AdjointMPSTensor{R₁}; trunc = notrunc(), kwargs...)

Left canonicalize a on-site MPS tensor.

# Return
-`::LeftIsometricTensor{R₁}`: left isometric tensor
-`::MPSTensor{2}`: bond tensor

If `trunc = notrunc()`, use `TensorKit.leftorth`, otherwise, use `TensorKit.tsvd`.
Propagate `kwargs` to the TensorKit functions.
"""
function leftorth(A::MPSTensor{R₁}; trunc=notrunc(), kwargs...) where R₁
    if trunc == notrunc()
        Q, R = leftorth(A.A, (Tuple(1:R₁-1), (R₁,)); kwargs...)
        if R₁ == 2 || R₁ == 3
            return LeftIsometricTensor(Q; check=false), MPSTensor(R), BondInfo(Q, :R)
        else
            return LeftIsometricTensor(permute(Q, (1, 2), Tuple(3:R₁)); check=false), MPSTensor(R), BondInfo(Q, :R)
        end
    else
        u, s, vd, info = tsvd(A, (Tuple(1:R₁-1), (R₁,)); trunc=trunc, kwargs...)
        if R₁ == 2 || R₁ == 3
            return LeftIsometricTensor(u; check=false), MPSTensor(s * vd), info
        else
            return LeftIsometricTensor(permute(u, (1, 2), Tuple(3:R₁)); check=false), MPSTensor(s * vd), info
        end
    end
end
function leftorth(A::AdjointMPSTensor{R₁}; trunc=notrunc(), kwargs...) where R₁
    if trunc == notrunc()
        if R₁ == 2
            Q, R = leftorth(A.A, ((2,), (1,)); kwargs...)
            return AdjointLeftIsometricTensor(permute(Q, (2,), (1,)); check=false), AdjointMPSTensor(permute(R, (2,), (1,))), BondInfo(Q, :R)
        else
            Q, R = leftorth(A.A, (Tuple(vcat(1:R₁-3, [R₁-1, R₁])), (R₁-2,)); kwargs...)
            return AdjointLeftIsometricTensor(permute(Q, Tuple(vcat(1:R₁-3, [R₁,])), (R₁-2, R₁-1)); check=false), AdjointMPSTensor(permute(R, (2,), (1,))), BondInfo(Q, :R)
        end
    else
        if R₁ == 2
            u, s, vd, info = tsvd(A, ((2,), (1,)); trunc=trunc, kwargs...)
            return AdjointLeftIsometricTensor(permute(u, (2,), (1,)); check=false), AdjointMPSTensor(permute(s * vd, (2,), (1,))), info
        else
            u, s, vd, info =  tsvd(A, (Tuple(vcat(1:R₁-3, [R₁-1, R₁])), (R₁-2,)); trunc=trunc, kwargs...)
            return AdjointLeftIsometricTensor(permute(u, Tuple(vcat(1:R₁-3, [R₁,])), (R₁-2, R₁-1)); check=false), AdjointMPSTensor(permute(s * vd, (2,), (1,))), info
        end
    end
end

"""
    rightorth{R₂}(A::MPSTensor; trunc = notrunc(), kwargs...)
    rightorth{R₂}(A::AdjointMPSTensor; trunc = notrunc(), kwargs...)

Right canonicalize a on-site MPS tensor.

# Return
-`::MPSTensor{2}`: bond tensor
-`::RightIsometricTensor{R₂}`: right isometric tensor

If `trunc = notrunc()`, use `TensorKit.rightorth`, otherwise, use `TensorKit.tsvd`.
Propagate `kwargs` to the TensorKit functions.
"""
function rightorth(A::MPSTensor{R₂}; trunc=notrunc(), kwargs...) where R₂
    if trunc == notrunc()
        L, Q = rightorth(A.A, ((1,), Tuple(2:R₂)); kwargs...)
        if R₂ == 2
            return MPSTensor(L), RightIsometricTensor(Q; check=false), BondInfo(Q, :L)
        else
            return MPSTensor(L), RightIsometricTensor(permute(Q, (1, 2), Tuple(3:R₂)); check=false), BondInfo(Q, :L)
        end
    else
        u, s, vd, info = tsvd(A, ((1,), Tuple(2:R₂)); trunc=trunc, kwargs...)
        if R₂ == 2
            return MPSTensor(u * s), RightIsometricTensor(vd; check=false), info
        else
            return MPSTensor(u * s), RightIsometricTensor(permute(vd, (1, 2), Tuple(3:R₂)); check=false), info
        end
    end
end
function rightorth(A::AdjointMPSTensor{R₂}; trunc=notrunc(), kwargs...) where R₂
    if trunc == notrunc()
        if R₂ == 2
            L, Q = rightorth(A.A, ((2,), (1, )); kwargs...)
            return MPSTensor(permute(L, (2,), (1,))), AdjointRightIsometricTensor(permute(Q, (2,), (1,)); check=false), BondInfo(Q, :L)
        else
            L, Q = rightorth(A.A, ((R₂-1,), Tuple(vcat((1:R₂-2), [R₂,]))); kwargs...)
            return MPSTensor(permute(L, (2,), (1,))), AdjointRightIsometricTensor(permute(Q, Tuple(2:R₂-1), (1, R₂)); check=false), BondInfo(Q, :L)
        end
    else
        if R₂ == 2
            u, s, vd, info = tsvd(A, ((2,), (1, )); trunc=trunc, kwargs...)
            return MPSTensor(permute(u * s, (2,), (1,))), AdjointRightIsometricTensor(permute(vd, (2,), (1,)); check=false), info
        else
            u, s, vd, info = tsvd(A, ((R₂-1,), Tuple(vcat((1:R₂-2), [R₂,]))); trunc=trunc, kwargs...)
            return MPSTensor(permute(u * s, (2,), (1,))), AdjointRightIsometricTensor(permute(vd, Tuple(2:R₂-1), (1, R₂)); check=false), info
        end
    end
end
