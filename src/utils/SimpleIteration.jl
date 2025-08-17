"""
    struct SimpleIteration
        tol::Union{Float64, Vector{Float64}}
        maxiter::Int64
    end

Simple iteration struct storing parameters

SimpleIteration is used to solve fixed point equation by iteration. See more for
    Base.iterate(f::Function, v₀::V; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter)
    Base.iterate(f::AbstractVector{<:Function}, v₀::AbstractVector{<:V}; tol::Vector{Float64}=repeat([Defaults.tol,], length(v₀)), maxiter::Int64=Defaults.maxiter)
    Base.iterate(f::AbstractVector{<:Function}, v₀::AbstractVector{<:V}, direction::Bool; tol::Vector{Float64}=repeat([Defaults.tol,], length(v₀)), maxiter::Int64=Defaults.maxiter)

# Constructor
SimpleIteration(; tol::Union{Float64, Vector{Float64}}=Defaults.tol, maxiter::Int64=Defaults.maxiter)
"""
struct SimpleIteration
    tol::Union{Float64, Vector{Float64}}
    maxiter::Int64

    function SimpleIteration(; tol::Union{Float64, Vector{Float64}}=Defaults.tol, maxiter::Int64=Defaults.maxiter)
        return new(tol, maxiter)
    end
end

const EigenAlgorithm = Union{SimpleIteration, KrylovKit.KrylovAlgorithm}
const GradientAlgorithm = Union{SimpleIteration}
