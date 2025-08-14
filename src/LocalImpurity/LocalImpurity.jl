"""
    abstract type AbstractLocalImpurity{W, L}

Abstract type of local impurity tensors, which correspond to local physical quantity.
"""
abstract type AbstractLocalImpurity{W, L} end

"""
    struct LocalImpurity{W, L}
        A::AbstractMatrix{MPOTensor}
    end

# Constructors
    LocalImpurity{W, L}(A::AbstractMatrix{<:MPOTensor})
    LocalImpurity(A::AbstractMatrix{<:MPOTensor})
"""
struct LocalImpurity{W, L} <: AbstractLocalImpurity{W, L}
    A::AbstractMatrix{MPOTensor}
    function LocalImpurity(A::AbstractMatrix{<:MPOTensor})
        W, L = size(A)
        return new{W, L}(A)
    end
end

convert(::Type{<:LocalImpurity}, A::AbstractMatrix{MPOTensor}) = LocalImpurity(A)

function show(io::IO, obj::LocalImpurity{W, L}) where {W, L}

    memory = obj |> Base.summarysize |> Base.format_bytes
    println(io, "$(typeof(obj)): total memory = $memory")
    # bond dimenson
    for w in 1:W
        println("line $w:")
        lsi = ceil(Int64, log10(L)) # length of si to be printed
        for si in 1:L
            D, DD = dim(obj[w, si], 1)
            print(io, "  Bond ", lpad(si-1, lsi), "->", lpad(si, lsi), ": ")
            println(io, "  $(sum(D)) -> $(sum(DD))")
        end
    end
end
