"""
    function Base.iterate(f::Function, v₀::V; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter) where V
    function Base.iterate(f::Function, v₀::V, alg::SimpleIteration) where V

# Arguments
`f::Function`: automorphism from V to V, i.e., for any v₀∈V, f(v₀)∈V.
`v₀::V`: initial value
`tol::Float64=Defaults.tol`: tolerance for iteration, i.e., norm(f(v)-v)
`maxiter::Int64=Defaults.maxiter`: iteration step limit

# Return
`v::V`: value after iteration
`info::SimpleIteratorInfo`: iteration information

# Functions to be defined
- Base.:(-)(::V, ::V) -> ::V
- norm(::V) -> Float64
"""
function Base.iterate(f::Function, v₀::V; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter) where V

    # information of simple iteration
    converged = false
    normres = Float64[]
    numiter = 0
    numops = 0

    v0 = deepcopy(v₀)

    while numiter < maxiter

        v = f(v0)
        numops += 1
        numiter += 1

        a = norm(v - v0)
        push!(normres, a)

        if a < tol
            converged = true
            break
        end

        v0 = deepcopy(v)
    end

    return v, SimpleIterationInfo(converged, normres, numiter, numops)
end
function Base.iterate(f::Function, v₀::V, alg::SimpleIteration) where V
    return iterate(f, v₀; tol=alg.tol, maxiter=alg.maxiter)
end
