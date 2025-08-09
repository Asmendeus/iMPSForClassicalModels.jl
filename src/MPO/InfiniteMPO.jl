"""
mutable struct InfiniteMPO{L, T} <: DenseInfiniteMPO{L}
    const A::AbstractVector{MPOTensor}
end

Concrete type of iMPO, where `L` is the cell size, `T == Float64` or `ComplexF64` is the number type of local tensors.

Graphic presentation:

           |            |
    ... -- A1 -- ... -- AL -- ...
           |            |

# Constructors
    InfiniteMPO{L, T}(::AbstractVector{<:MPOTensor})
    InfiniteMPO(::AbstractVector{<:MPOTensor})
    InfiniteMPO{L, T}()
    InfiniteMPO(L, T=Float64)

Note:
    `InfiniteMPO` and its adjoint are similar to `InfiniteMPS`, used only as a variational boundary, and the partition function corresponds to `SparseMPO`.
    There are two kinds of fixed point equations for a MPO boundary:

        1.
               |    |    |    |                                          |    |    |    |
            -- A -- A -- A -- A --          |    |    |    |          —— O —— O —— O —— O ——          |    |    |    |
               |    |    |    |     =  λ -- A -- A -- A -- A --  and     |    |    |    |     =  λ -- A -- A -- A -- A --
            —— O —— O —— O —— O ——          |    |    |    |          -- A -- A -- A -- A --          |    |    |    |
               |    |    |    |                                          |    |    |    |

        This case means that we generalize 2D classical system with periodic boundary on one direction and one open boundary on the other to the thermodynamic limit.

        2.
               |    |    |    |                                          |    |    |    |
            -- A -- A -- A -- A --          |    |    |    |          —— O —— O —— O —— O ——          |    |    |    |
               |    |    |    |     =  λ -- A -- A -- A -- A --  and     |    |    |    |     =  λ -- B -- B -- B -- B --
            —— O —— O —— O —— O ——          |    |    |    |          -- B -- B -- B -- B --          |    |    |    |
               |    |    |    |                                          |    |    |    |

        This case means that we map the 2D classical system to a 1D quantum system, whose ground state is a mixed state, if the following rank-1 decompositions can not be strictly implemented

                           |                            |
               |           E                |        -- F --
            -- A --  =     ⊗            -- B --  =     ⊗
               |        -- D --             |           G
                           |                            |

    I have no idea if these two cases will have the same results as the case with MPS as the boundary, but I think the same may be the norm.
    Anyway, this package support MPO as a variational boundary for the 2D system.
    If you have any interesting ideas or published results on these two or other cases where variational boundary should be a MPO, welcome to discuss with me. I'm interested in this.
    Incidentally, promoting the 2D classical model with periodic boundaries in both directions to the thermodynamic limit implies promoting the environmen tensor, whose computational complexity is approximately the square of the existing algorithm.
    This package currently has no implementation for this case.
"""
mutable struct InfiniteMPO{L, T} <: DenseInfiniteMPO{L}
    const A::AbstractVector{MPOTensor}

    function InfiniteMPO{L, T}(A::AbstractVector{<:MPOTensor}) where {L, T}
        length(A) == L || throw(ArgumentError("Mismatched lengths: $L ≠ $(length(A))"))
        T ∈ [Float64, ComplexF64] || throw(ArgumentError("Unsupported data types $T"))
        return new{L, T}(A)
    end
    function InfiniteMPO(A::AbstractVector{<:MPOTensor})
        L = length(A)
        T = mapreduce(eltype, promote_type, A)
        return new{L, T}(A)
    end

    function InfiniteMPO{L, T}() where {L, T}
        A = Vector{MPOTensor}(undef, L)
        return new{L, T}(A)
    end
    InfiniteMPO(L::Int, T::Type{<:Union{Float64, ComplexF64}}=Float64)=InfiniteMPO{L, T}()
end

const iMPO = InfiniteMPO
