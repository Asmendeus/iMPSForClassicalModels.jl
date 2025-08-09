"""
    mutable struct IsometricEnvironment{N, R} <: AbstractEnvironment{N}
        const A::LocalTensor{R}
        const O::AbstractVector{MPOTensor}
        const B::AdjointLocalTensor{R}
    end

Wrapper type for left environment tensor's generating environment.

Graphic presentation:

    -- A  --    -- A  --
       |           |
    —— O1 ——    —— O1 ——
       |           |
       ⋮            ⋮
       |           |
    —— OW ——    —— OW ——
       |           |
    -- B  --    -- B  --

# Constructors
    IsometricEnvironment{N, R}(A::LocalTensor{R}, O::AbstractVector{<:MPOTensor}, B::AdjointLocalTensor{R})
    IsometricEnvironment(A::LocalTensor{R}, O::AbstractVector{<:MPOTensor}, B::AdjointLocalTensor{R})
    IsometricEnvironment(A::LocalTensor{R}, O::MPOTensor, B::AdjointLocalTensor{R})
"""
mutable struct IsometricEnvironment{N, R} <: AbstractEnvironment{N}
    A::LocalTensor{R}
    const O::AbstractVector{MPOTensor}
    B::AdjointLocalTensor{R}

    function IsometricEnvironment{N, R}(A::LocalTensor{R}, O::AbstractVector{<:MPOTensor}, B::AdjointLocalTensor{R}) where {N, R}
        N == length(O) + 2 || throw(ArgumentError("Mismatched widths: $N ≠ $(length(O)+2)"))
        R == 2 && throw(ArgumentError("Not support bond tensors in forming left environment"))

        (isLeftIsometric(A) || isRightIsometric(A)) || throw(ArgumentError("The local tensor is not left or right orthogonal"))
        (isLeftIsometric(B) || isRightIsometric(B)) || throw(ArgumentError("The adjoint local tensor is not left or right orthogonal"))

        pspace_Adn = space(A, 2)
        pspace_Oup = space(O[1], 3)
        pspace_Adn == pspace_Oup || throw(SpaceMismatch("Mismatched physical spaces of local tensor and MPO tensor: $(pspace_Adn) ≠ $(pspace_Oup))"))

        for w in 1:N-3
            pspace_Odn_w = space(O[w], 2)
            pspace_Oup_w = space(O[w+1], 3)
            pspace_Odn_w == pspace_Oup_w || throw(SpaceMismatch("Mismatched physical spaces of $w and $(w+1) layer MPO tensors: $(pspace_Odn_w) ≠ $(pspace_Oup_w))"))
        end

        pspace_Bup = space(B, R)
        pspace_Odn = space(O[end], 2)
        pspace_Bup == pspace_Odn || throw(SpaceMismatch("Mismatched physical spaces of MPO tensor and adjoint local tensor: $(pspace_Bup) ≠ $(pspace_Odn))"))

        if R > 3
            pspace_Aup = [space(A, i) for i in 3:R-1]
            pspace_Bdn = [space(B, i) for i in 1:R-3]
            pspace_Aup == pspace_Bdn || throw(SpaceMismatch("Mismatched physical spaces of local tensor and adjoint local tensor: $(pspace_Bup) ≠ $(pspace_Odn))"))
        end

        return new{N, R}(A, O, B)
    end
    IsometricEnvironment(A::LocalTensor{R}, O::AbstractVector{<:MPOTensor}, B::AdjointLocalTensor{R}) where R = IsometricEnvironment{length(O)+2, R}(A, O, B)
    IsometricEnvironment(A::LocalTensor{R}, O::MPOTensor, B::AdjointLocalTensor{R}) where R = IsometricEnvironment{3, R}(A, [O,], B)
end

const IsometricMPSEnvironment{N} = IsometricEnvironment{N, 3}
const IsometricMPOEnvironment{N} = IsometricEnvironment{N, 4}
