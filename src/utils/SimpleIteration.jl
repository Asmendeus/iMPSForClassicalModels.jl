"""
    struct SimpleIteration
        tol::Float64
        maxiter::Int64
    end

SimpleIteration is solving fixed point equation by iteration. See more for
    Base.iterate(f::Function, v₀::V; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter)
    Base.iterate(f::Function, v₀::V, alg::SimpleIteration)

# Constructor
SimpleIteration(; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter)
"""
struct SimpleIteration
    tol::Float64
    maxiter::Int64

    function SimpleIteration(; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter)
        return new(tol, maxiter)
    end
end
