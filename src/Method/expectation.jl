function expecatation(mps1::InfiniteMPS{L, T}, mps2::AdjointInfiniteMPS{L, T}, mpo::InfiniteMPO{L, T}, M::LocalImpurity{L}) where {L, T}
end
function expecatation(mps1::InfiniteMPS{L, T}, mps2::AdjointInfiniteMPS{L, T}, mpo::MultirowInfiniteMPO{L, W, T}, M::MultirowLocalImpurity{L, W}) where{L, W, T}
end
function expecatation(mps1::InfiniteMPS{L, T}, mps2::AdjointInfiniteMPS{L, T}, mpo::InfiniteMPO{L, T}, M1::LocalImpurity{L}, M2::LocalImpurity{L}, r::Int64) where {L, T}
end
function expecatation(mps1::InfiniteMPS{L, T}, mps2::AdjointInfiniteMPS{L, T}, mpo::MultirowInfiniteMPO{L, W, T}, M1::MultirowLocalImpurity{L, W}, M2::MultirowLocalImpurity{L, W}, r::Int64) where{L, W, T}
end