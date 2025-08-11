"""
    struct SimpleIteration
        tol::Union{Float64, Vector{Float64}}
        maxiter::Int64
    end

SimpleIteration is solving fixed point equation by iteration. See more for
    Base.iterate(f::Function, v₀::V; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter)
    Base.iterate(f::Function, v₀::Vector{<:V}; tol::Vector{Float64}=repeat([Defaults.tol,], length(v₀)), maxiter::Int64=Defaults.maxiter)
    Base.iterate(f::Function, v₀::V, alg::SimpleIteration)

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
