"""
    Base.iterate(f::Function, v₀::V; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter)
    Base.iterate(f::Function, v₀::V, alg::SimpleIterator)

A simple iterative procedure for solving fixed point equation
    f(v) = λ * v,
where `f` is a self-map, i.e., `f(::V) -> ::V`. `λ` is normalization coefficient of fixed point `v`.

# Arguments
`f::Function`: self-map on `V`, i.e., for any v∈V, f(v)∈V
`v₀::V`: initial argument
`alg::SimpleIterator`: wrapper for keyword arguments

# Keyword Arguments
`tol::Float64=Defaults.tol`: tolerance for iteration
`maxiter::Int64=Defaults.maxiter`: iteration step limit

# Return
`λ::N`: normalization coefficient of `f(v)`
`v::V`: argument after iteration
`info::SimpleIteratorInfo`: iteration information

# Functions to be defined
`lambda(::V) -> ::N`
`division(::V, ::N) -> ::V`
`minus(::V, ::V) -> ::V`
`norm_max(::V) -> Float64`
"""
function Base.iterate(f::Function, v₀::V; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter) where V

    converged = false
    normres = Float64[]
    numiter = 0
    numops = 0

    λ = lambda(v₀)
    v₀ = division(v₀, λ)

    while numiter < maxiter

        v = f(v₀)

        λ = lambda(v)
        v = division(v, λ)

        numops += 1
        numiter += 1

        nr = norm_max(minus(v, v₀))
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

    return λ, v₀, SimpleIteratorInfo(converged, normres, numiter, numops)
end
function Base.iterate(f::Function, v₀::V, alg::SimpleIterator) where V
    return iterate(f, v₀; tol=alg.tol, maxiter=alg.maxiter)
end

"""
    Base.iterate(f::AbstractVector{<:Function}, v₀::AbstractVector{<:V}; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter)
    Base.iterate(f::AbstractVector{<:Function}, v₀::AbstractVector{<:V}, alg::SimpleIterator)

A simple iterative procedure for solving fixed point equations

    v[1] = f[1](v[1]) / λ[1]
    v[2] = f[2](v[2]) / λ[2]
    v[3] = f[3](v[3]) / λ[3]
            ⋮
    v[end] = f[end](v[end]) / λ[end]

where `f` are self-maps, i.e., `f[i](::V) -> ::V`. `λ` are normalization coefficients of fixed points `v`.

# Arguments
`f::AbstractVector{<:Function}`: self-maps on `V`, i.e., for any v∈V, f[i](v)∈V
`v₀::AbstractVector{<:V}`: initial arguments
`alg::SimpleIterator`: wrapper for keyword arguments

# Keyword Arguments
`tol::Float64=Defaults.tol`: tolerance for iteration
`maxiter::Int64=Defaults.maxiter`: iteration step limit

# Return
`v::V`: arguments after iteration
`info::SimpleIteratorInfo`: iteration information

# Functions to be defined
`lambda(::AbstractVector{<:V}) -> ::AbstractVector{<:N}`
`division(::AbstractVector{<:V}, ::AbstractVector{<:N}) -> ::AbstractVector{<:V}`
`minus(::AbstractVector{<:V}, ::AbstractVector{<:V}) -> ::AbstractVector{<:V}`
`norm_max(::AbstractVector{<:V}) -> Float64`
"""
function Base.iterate(f::AbstractVector{<:Function}, v₀::AbstractVector{<:V}; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter) where V

    (L = length(f)) == length(v₀) || throw(ArgumentError("Mismatched lengths of `f` and `v₀`: ($(length(f)), $(length(v₀)))"))

    converged = false
    normres = Float64[]
    numiter = 0
    numops = 0

    λ = lambda(v₀)
    v₀ = division(v₀, λ)

    while numiter < maxiter

        v = map(l->f[l](v₀[l]), 1:L)

        λ = lambda(v)
        v = division(v, λ)

        numops += 1
        numiter += 1

        nr = norm_max(minus(v, v₀))
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

    return λ, v₀, SimpleIteratorInfo(converged, normres, numiter, numops)
end
function Base.iterate(f::AbstractVector{<:Function}, v₀::AbstractVector{<:V}, alg::SimpleIterator) where V
    return iterate(f, v₀; tol=alg.tol, maxiter=alg.maxiter)
end

"""
    Base.iterate(f::AbstractVector{<:Function}, v₀::AbstractVector{<:V}, direction::Bool; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter)
    Base.iterate(f::AbstractVector{<:Function}, v₀::AbstractVector{<:V}, direction::Bool, alg::SimpleIterator)

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

where `f` are self-maps, i.e., `f[i](::V) -> ::V`. `λ` are normalization coefficients of fixed points `v`.

# Arguments
`f::AbstractVector{<:Function}`: self-maps on `V`, i.e., for any v∈V, f[i](v)∈V
`v₀::AbstractVector{<:V}`: initial arguments
`direction::Bool`: loop direction - true is positive loop and false is reverse loop
`alg::SimpleIterator`: wrapper for keyword arguments

# Keyword Arguments
`tol::Float64=Defaults.tol`: tolerance for iteration
`maxiter::Int64=Defaults.maxiter`: iteration step limit

# Return
`v::V`: arguments after iteration
`info::SimpleIteratorInfo`: iteration information

# Functions to be defined
`lambda(::AbstractVector{<:V}) -> ::AbstractVector{<:N}`
`division(::AbstractVector{<:V}, ::AbstractVector{<:N}) -> ::AbstractVector{<:V}`
`minus(::AbstractVector{<:V}, ::AbstractVector{<:V}) -> ::AbstractVector{<:V}`
`norm_max(::AbstractVector{<:V}) -> Float64`
"""
function Base.iterate(f::AbstractVector{<:Function}, v₀::AbstractVector{<:V}, direction::Bool; tol::Float64=Defaults.tol, maxiter::Int64=Defaults.maxiter) where V

    (L = length(f)) == length(v₀) || throw(ArgumentError("Mismatched lengths of `f` and `v₀`: ($(length(f)), $(length(v₀)))"))

    converged = false
    normres = Float64[]
    numiter = 0
    numops = 0

    λ = lambda(v₀)
    v₀ = division(v₀, λ)
    v = deepcopy(v₀)

    while numiter < maxiter

        if direction
            for l in 1:L
                l₊ = mod(l, L)+1

                v[l₊] = f[l](v₀[l])

                λ[l] = lambda(v[l₊])
                v[l₊] /= λ[l]
            end
        else
            for l in L:-1:1
                l₋ = mod(l-2, L) + 1

                v[l₋] = f[l](v₀[l])

                λ[l] = lambda(v[l₋])
                v[l₋] /= λ[l]
            end
        end

        numops += L
        numiter += 1

        nr = norm_max((minus(v, v₀)))
        push!(normres, nr)

        v₀ = deepcopy(v)

        if nr < tol
            converged = true
            break
        end
    end

    if !converged
        @warn "Simple iterations fail to converge: ϵ = $(normres[end]) .< $tol not all true"
    end

    return λ, v₀, SimpleIteratorInfo(converged, normres, numiter, numops)
end
function Base.iterate(f::AbstractVector{<:Function}, v₀::AbstractVector{<:V}, direction::Bool, alg::SimpleIterator) where V
    return iterate(f, v₀, direction; tol=alg.tol, maxiter=alg.maxiter)
end
