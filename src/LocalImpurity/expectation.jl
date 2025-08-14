function expecatation(ψ::DenseInfiniteMPS{L}, H::SparseMPO{W, L}, ψ′::AdjointInfiniteMPS{L}, M::LocalImpurity{W, L};
            FL₀::AbstractVector{<:LeftEnvironmentTensor{W+2}},
            FR₀::AbstractVector{<:RightEnvironmentTensor{W+2}},
            alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()) where {W, L}
end

function expecatation(ψ::DenseInfiniteMPS{L}, H::SparseMPO{W, L}, ψ′::AdjointInfiniteMPS{L},
            M1::LocalImpurity{W, L}, M2::LocalImpurity{W, L}, r::Int64;
            FL₀::AbstractVector{<:LeftEnvironmentTensor{W+2}},
            FR₀::AbstractVector{<:RightEnvironmentTensor{W+2}},
            alg::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}=Arnoldi()) where {W, L}
end