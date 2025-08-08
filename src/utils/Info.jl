"""
     struct BondInfo
          D::Int64
          DD::Int64
          TrunErr::Float64
          SE::Float64
     end

Type for storing the information of a bond.

# Constructors
     BondInfo(s::AbstractTensorMap, ϵ::Float64 = 0.0)

Outer constructor via giving the `s` tensor and `ϵ` form `tsvd`.

     BondInfo(A::AbstractTensorMap, direction::Symbol)
     BondInfo(A::MPSTensor, direction::Symbol)

Outer constructor via giving a tensor `A` and `direction = :L` or `:R`. We cannot get truncation error and singular values hence `TrunErr` and `SE` are set to `0.0` and `NaN`, respectively.
"""
struct BondInfo
     D::Int64
     DD::Int64
     TrunErr::Float64
     SE::Float64
     BondInfo(D::Int64, DD::Int64, TrunErr::Float64, SE::Float64) = new(D, DD, TrunErr, SE)
end

function BondInfo(s::AbstractTensorMap{<:Union{Float64, ComplexF64} ,T}, ϵ::Float64=0.0) where T <: GradedSpace
     D = DD = 0
     Norm2 = SE = 0.0
     for (c, b) in blocks(s)
          λ = diag(b)
          D += length(λ)
          DD += length(λ) * dim(c)
          Norm2 += norm(λ)^2 * dim(c)
          SE += mapreduce(x -> x == 0 ? 0 : x^2 * log(x), +, λ) * dim(c)
     end
     SE = -2SE / Norm2 + log(Norm2)
     return BondInfo(D, DD, ϵ, SE)
end
function BondInfo(s::AbstractTensorMap{<:Union{Float64, ComplexF64}, T}, ϵ::Float64=0.0) where T <: Union{CartesianSpace, ComplexSpace}
    D = DD = 0
    Norm2 = SE = 0.0

    λ = data(s)[1]
    D += length(λ)
    DD += length(λ)
    Norm2 += norm(λ)^2
    SE += mapreduce(x -> x == 0 ? 0 : x^2 * log(x), +, λ; init = 0.0)

    SE = -2SE / Norm2 + log(Norm2)
    return BondInfo(D, DD, ϵ, SE)
end
function BondInfo(A::AbstractTensorMap, direction::Symbol)
     @assert direction in (:L, :R)
     idx = direction == :L ? 1 : numind(A)
     return BondInfo(dim(A, idx)..., 0.0, NaN)
end

function show(io::IO, info::BondInfo)
     print(io, "BondInfo(D = $(info.D) => $(info.DD), TrunErr2 = $(info.TrunErr^2), SE = $(info.SE))")
end

function merge(info1::BondInfo, info2::BondInfo)

     if isnan(info1.SE)
          SE = info2.SE
     elseif isnan(info2.SE)
          SE = info1.SE
     else
          SE = max(info1.SE, info2.SE)
     end

     return BondInfo(max(info1.D, info2.D),
          max(info1.DD, info2.DD),
          max(info1.TrunErr, info2.TrunErr),
          SE)
end
merge(info1::BondInfo, info2::BondInfo, args...) = merge(merge(info1, info2), args...)
merge(v::AbstractVector{BondInfo}) = reduce(merge, v)
