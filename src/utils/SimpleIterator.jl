"""
    struct SimpleIterator
        tol::Float64
        maxiter::Int64
    end

Simple iteration struct storing parameters

SimpleIterator is used to solve fixed point equation by iteration. See more for
    Base.iterate(f::Function, v₀::V; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter)
    Base.iterate(f::AbstractVector{<:Function}, v₀::AbstractVector{<:V}; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter)
    Base.iterate(f::AbstractVector{<:Function}, v₀::AbstractVector{<:V}, direction::Bool; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter)

# Constructor
    SimpleIterator(; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter)
"""
struct SimpleIterator
    tol::Float64
    maxiter::Int64

    function SimpleIterator(; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter)
        return new(tol, maxiter)
    end
end

const EigenAlgorithm = Union{SimpleIterator, KrylovKit.KrylovAlgorithm}
const GradientAlgorithm = Union{SimpleIterator, OptimKit.OptimizationAlgorithm}
