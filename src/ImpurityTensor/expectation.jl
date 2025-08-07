function expecatation(mps1::InfiniteMPS{L, T}, mps2::AdjointInfiniteMPS{L, T}, mpo::InfiniteMPO{L, T}, M::ImpurityTensor{L}) where {L, T}
end
function expecatation(mps1::InfiniteMPS{L, T}, mps2::AdjointInfiniteMPS{L, T}, mpo::MultirowInfiniteMPO{L, W, T}, M::MultirowImpurityTensor{L, W}) where{L, W, T}
end
function expecatation(mps1::InfiniteMPS{L, T}, mps2::AdjointInfiniteMPS{L, T}, mpo::InfiniteMPO{L, T}, M1::ImpurityTensor{L}, M2::ImpurityTensor{L}, r::Int64) where {L, T}
end
function expecatation(mps1::InfiniteMPS{L, T}, mps2::AdjointInfiniteMPS{L, T}, mpo::MultirowInfiniteMPO{L, W, T}, M1::MultirowImpurityTensor{L, W}, M2::MultirowImpurityTensor{L, W}, r::Int64) where{L, W, T}
end