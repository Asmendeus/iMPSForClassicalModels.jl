"""
    mutable struct ChannelEnvironment{N, L, R} <: AbstractEnvironment{N}
        const A::AbstractVector{<:LocalTensor{R}}
        const O::AbstractMatrix{MPOTensor}
        const B::AbstractVector{<:AdjointLocalTensor{R}}
    end

Wrapper type for left/right environment tensor's generating environment.

Graphic presentation:
                                        (a)   (bonds `a` are connected)
                                         |
    -- A[1] ---- ... -- A[L] ----     -- A[1] ---- ... -- A[L] ----
       |                |                |                |
    —— O[1,1] —— ... —— O[1,L] ——     —— O[1,1] —— ... —— O[1,L] ——
       |                |                |                |
       ⋮                 ⋮                ⋮                 ⋮
       |                |                |                |
    —— O[W,1] —— ... —— O[W,L] ——     —— O[W,L] —— ... —— O[W,L] ——  W = N - 2
       |                |                |                |
    -- B[1] ---- ... -- B[L] ----     -- B[1] ---- ... -- B[L] ----
                                         |
                                        (a)

# Constructors
    ChannelEnvironment{N, L, R}(A::AbstractVector{<:LocalTensor{R}}, O::AbstractMatrix{<:MPOTensor}, B::AbstractVector{<:AdjointLocalTensor{R}})
    ChannelEnvironment(A::LocalTensor{R}, O::AbstractVector{<:MPOTensor}, B::AdjointLocalTensor{R})
    ChannelEnvironment(A::LocalTensor{R}, O::MPOTensor, B::AdjointLocalTensor{R})
"""
mutable struct ChannelEnvironment{N, L, R} <: AbstractEnvironment{N}
    const A::AbstractVector{<:LocalTensor{R}}
    const O::AbstractMatrix{MPOTensor}
    const B::AbstractVector{<:AdjointLocalTensor{R}}

    function ChannelEnvironment{N, L, R}(A::AbstractVector{<:LocalTensor{R}}, O::AbstractMatrix{<:MPOTensor}, B::AbstractVector{<:AdjointLocalTensor{R}}) where {N, L, R}
        N == size(O, 1) + 2 || throw(ArgumentError("Mismatched widths: $N ≠ $(size(O, 1)+2)"))
        L == length(A) == size(O, 2) == length(B) || throw(ArgumentError("Mismatched lengths: ($L, $(length(A)), $(size(O, 2)), $(length(B)))"))
        R == 2 && throw(ArgumentError("Not support bond tensors in forming channel environment"))

        for l in 1:L
            pspace_Aup_l = [domain(A[l], i) for i in 1:R-3]
            pspace_Bdn_l = [codomain(B[l], i) for i in 1:R-3]
            pspace_Aup_l == pspace_Bdn_l || throw(SpaceMismatch("Mismatched physical space between $(l)-th local tensor and adjoin local tensor: $(pspace_Aup_l) ≠ $(pspace_Bdn_l))"))

            pspace_Adn_l = codomain(A[l], 2)
            pspace_O1up_l = domain(O[1, l], 1)
            pspace_Adn_l == pspace_O1up_l || throw(SpaceMismatch("Mismatched physical space between $(l)-th local tensor and mpo tensor: $(pspace_Adn_l) ≠ $(pspace_O1up_l)"))

            pspace_Bup_l = domain(B[l], 2)
            pspace_OWdn_l = codomain(O[N-2, l], 2)
            pspace_Bup_l == pspace_OWdn_l || throw(SpaceMismatch("Mismatched physical space between $(l)-th mpo tensor and adjoint local tensor: $(pspace_Bup_l) ≠ $(pspace_OWdn_l)"))

            for w in 1:N-3
                pspace_Owdn_l = codomain(O[w, l], 2)
                pspace_Owup_l = domain(O[w+1, l], 1)
                pspace_Owdn_l == pspace_Owup_l || throw(SpaceMismatch("Mismatched physical space between ($w, $l)-th and ($(w+1), $l)-th mpo tensors: $(pspace_Owdn_l) ≠ $(pspace_Owup_l)"))
            end

            aspace_Art_l = domain(A[l], R-2)
            aspace_Alt_l = codomain(A[mod(l, L)+1], 1)
            aspace_Art_l == aspace_Alt_l || throw(SpaceMismatch("Mismatched virtual space between $(l)-th and $(mod(l, L)+1)-th local tensors: $(aspace_Art_l) ≠ $(aspace_Alt_l))"))

            aspace_Brt_l = codomain(B[l], R-2)
            aspace_Blt_l = domain(B[mod(l, L)+1], 1)
            aspace_Brt_l == aspace_Blt_l || throw(SpaceMismatch("Mismatched virtual space between $(l)-th and $(mod(l, L)+1)-th adjoint local tensors: $(aspace_Brt_l) ≠ $(aspace_Blt_l))"))

            for w in 1:N-2
                aspace_Owrt_l = domain(O[w, l], 2)
                aspace_Owlt_l = codomain(O[w, mod(l, L)+1], 1)
                aspace_Owrt_l == aspace_Owlt_l || throw(SpaceMismatch("Mismatched auxilliary space between ($w, $l)-th and ($w, $(mod(l, L)+1))-th mpo tensors: $(aspace_Owrt_l) ≠ $(aspace_Owlt_l)"))
            end
        end

        return new{N, L, R}(A, O, B)
    end
    ChannelEnvironment(A::AbstractVector{<:LocalTensor{R}}, O::AbstractMatrix{<:MPOTensor}, B::AbstractVector{<:AdjointLocalTensor{R}}) where R = ChannelEnvironment{size(O, 1)+2, length(A), R}(A, O, B)
    ChannelEnvironment(A::LocalTensor{R}, O::AbstractVector{<:MPOTensor}, B::AdjointLocalTensor{R}) where R = ChannelEnvironment{length(O)+2, 1, R}([A,], reshape(O, :, 1), [B,])
    ChannelEnvironment(A::LocalTensor{R}, O::MPOTensor, B::AdjointLocalTensor{R}) where R = ChannelEnvironment{3, 1, R}([A,], [O;;], [B,])
end

const MPSChannelEnvironment{N, L} = ChannelEnvironment{N, L, 3}
const MPOChannelEnvironment{N, L} = ChannelEnvironment{N, L, 4}
