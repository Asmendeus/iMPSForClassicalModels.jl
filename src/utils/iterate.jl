"""
    Base.iterate(f::Function, x₀::Tuple{Number, V}; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter) where V
    Base.iterate(f::Function, v₀::V; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter) where V
    Base.iterate(f::Function, v₀::V, alg::SimpleIteration) where V

A simple iterative procedure for solving fixed point equation
    v = f(v) / λ,
where `λ` is normalization coefficient of fixed point with `norm(v) = 1` and `norm(f(v)) = λ`.

# Arguments
`f::Function`: automorphism from V to V, i.e., for any v∈V, f(v)∈V
`x₀::Tuple{Number, V}`: x₀ = (_, v₀)
`v₀::V`: initial value
`alg::SimpleIteration`: wrapper for keyword arguments

# Keyword Arguments
`tol::Float64=Defaults.tol`: tolerance for iteration, i.e., norm(f(v) / λ - v)
`maxiter::Int64=Defaults.maxiter`: iteration step limit

# Return
`v::V`: value after iteration
`info::SimpleIterationInfo`: iteration information

# Functions to be defined
- Base.:(-)(::V, ::V) -> ::V
- Base.:(/)(::V, ::Float64) -> ::V
- norm(::V) -> Float64
- sign_first_element(::V) -> Number
"""
function Base.iterate(f::Function, x₀::Tuple{Number, V}; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter) where V

    converged = false
    normres = Float64[]
    numiter = 0
    numops = 0

    λ = norm(x₀[2]) * sign_first_element(x₀[2])
    v₀ = x₀[2] / λ

    while numiter < maxiter

        v = f(v₀)

        λ = norm(v) * sign_first_element(v)
        v /= λ

        numops += 1
        numiter += 1

        nr = norm(v - v₀)
        push!(normres, nr)

        v₀ = deepcopy(v)

        if nr < tol
            converged = true
            break
        end
    end

    if !converged
        @warn "Simple iterations fail to converge: ϵ = $(normres[end]) > $tol"
    end

    return λ, v₀, SimpleIterationInfo(converged, normres, numiter, numops)
end
function Base.iterate(f::Function, v₀::V; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter) where V
    return iterate(f, (1.0, v₀); tol=tol, maxiter=maxiter)
end
function Base.iterate(f::Function, v₀::V, alg::SimpleIteration) where V
    return iterate(f, (1.0, v₀); tol=alg.tol, maxiter=alg.maxiter)
end

"""
    Base.iterate(f::Vector{<:Function}, x₀::Tuple{Vector{<:Number}, Vector{<:V}}, direction::Bool; tol::Vector{Float64}=repeat([Defaults.tol,], length(f)), maxiter::Int64=Defaults.maxiter) where V
    Base.iterate(f::Vector{<:Function}, v₀::Vector{<:V}, direction::Bool; tol::Vector{Float64}=repeat([Defaults.tol,], length(v₀)), maxiter::Int64=Defaults.maxiter) where V
    Base.iterate(f::Vector{<:Function}, v₀::Vector{<:V}, direction::Bool, alg::SimpleIteration) where V

A simple iterative procedure for solving loop fixed point equations

    if direction == true (positive loop)
        v[2] = f[1](v[1]) / λ[1]
        v[3] = f[2](v[2]) / λ[2]
                ⋮
        v[end] = f[end-1](v[end-1]) / λ[end-1]
        v[1] = f[end](v[end]) / λ[end]

    elseif direction == false (reverse loop)
        v[end-1] = f[end](v[end]) / λ[end]
        v[end-2] = f[end-1](v[end-1]) / λ[end-1]
                ⋮
        v[1] = f[2](v[2]) / λ[2]
        v[end] = f[1](v[1]) / λ[1]

where `λ` are normalization coefficients of fixed point with `norm(v[i]) = 1` and `norm(f[i](v[i])) = λ[i]`.

# Arguments
`f::Vector{<:Function}`: automorphisms from V to V, i.e., for any v∈V, f(v)∈V
`x₀::Tuple{Vector{<:Number}, Vector{<:V}}`: x₀ = (_, v₀)
`v₀::Vector{<:V}`: initial values
`direction::Bool`: loop direction - true is positive loop and false is reverse loop
`alg::SimpleIteration`: wrapper for keyword arguments

# Keyword Arguments
`tol::Vector{Float64}=repeat([Defaults.tol,], length(x₀[2]))`: tolerances for iteration, i.e., norm(f.(v) ./ λ - v)
`maxiter::Int64=Defaults.maxiter`: iteration step limit

# Return
`v::V`: values after iteration
`info::SimpleIterationInfo`: iteration information

# Functions to be defined
- Base.:(-)(::V, ::V) -> ::V
- Base.:(/)(::V, ::Float64) -> ::V
- norm(::V) -> Float64
- sign_first_element(::V) -> Number
"""
function Base.iterate(f::Vector{<:Function}, x₀::Tuple{Vector{<:Number}, Vector{<:V}}, direction::Bool; tol::Vector{Float64}=repeat([Defaults.tol,], length(f)), maxiter::Int64=Defaults.maxiter) where V

    (L = length(f)) == length(x₀[2]) == length(tol) || throw(ArgumentError("Mismatched lengths of `f`, `v₀` and `tol`: ($(length(f)), $(length(x₀[2])), $(length(tol)))"))

    converged = false
    normres = Vector{Float64}[]
    numiter = 0
    numops = 0

    λ = norm.(x₀[2]) .* sign_first_element.(x₀[2])
    v₀ = x₀[2] ./ λ
    v = deepcopy(v₀)

    while numiter < maxiter

        if direction
            for l in 1:L
                l₊ = mod(l, L)+1

                v[l₊] = f[l](v₀[l])

                λ[l] = norm(v[l₊]) * sign_first_element(v[l₊])
                v[l₊] /= λ[l]
            end
        else
            for l in L:-1:1
                l₋ = mod(l-2, L) + 1

                v[l₋] = f[l](v₀[l])

                λ[l] = norm(v[l₋]) * sign_first_element(v[l₋])
                v[l₋] /= λ[l]
            end
        end

        numops += L
        numiter += 1

        nr = norm.(v - v₀)
        push!(normres, nr)

        v₀ = deepcopy(v)

        if all(nr .< tol)
            converged = true
            break
        end
    end

    if !converged
        @warn "Simple iterations fail to converge: ϵ = $(normres[end]) .< $tol not all true"
    end

    return λ, v₀, SimpleIterationInfo(converged, normres, numiter, numops)
end
function Base.iterate(f::Vector{<:Function}, v₀::Vector{<:V}, direction::Bool; tol::Vector{Float64}=repeat([Defaults.tol,], length(v₀)), maxiter::Int64=Defaults.maxiter) where V
    return iterate(f, ([1.0,], v₀), direction; tol=tol, maxiter=maxiter)
end
function Base.iterate(f::Vector{<:Function}, v₀::Vector{<:V}, direction::Bool, alg::SimpleIteration) where V
    return iterate(f, ([1.0,], v₀), direction; tol=alg.tol, maxiter=alg.maxiter)
end

# Normalize eigenvectors
sign_first_element(A::Number) = sign(A)
sign_first_element(A::AbstractArray) = sign(A[1])
sign_first_element(A::AbstractTensorMap) = sign(A[1])
# sign_first_element(A::AbstractTensorWrapper) = sign(A.A[1])
