"""
    function isLeftIsometric(A::AbstractTensorMap; tol::Float64=1e-8)
    function isLeftIsometric(A::AbstractMPSTensor; tol::Float64=1e-8)

# Return
    `::Bool`: whether the 3-leg MPS tensor `A` is left orthogonal
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
    `::Bool`: whether the 3-leg MPS tensor `A` is right orthogonal
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
    function leftorth(A::MPSTensor)
"""
function leftorth(A::MPSTensor)
end

"""
    function rightorth(A::MPSTensor)
"""
function rightorth(A::MPSTensor)
end

