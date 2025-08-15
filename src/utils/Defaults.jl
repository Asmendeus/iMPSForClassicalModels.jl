module Defaults
import KrylovKit.Arnoldi
include("SimpleIteration.jl")

    datatype = ComplexF64

    dense = 0.1         # sparse matrix is matrix with dense degree smaller than it

    alg_eig = Arnoldi           # algorithm for solving maximum/minimum eigenvalue problem
    alg_grad = SimpleIteration  # algorithm for solving gradient problem

    tol_high = 1e-16        # high-accuracy iteration tolerance
    tol = 1e-14             # normal-accuracy iteration tolerance
    tol_low = 1e-12         # low-accuracy iteration tolerance

    maxiter = 200       # iteration step limit

    tol_norm = 1e-10    # norm(A - B): tolerance of equation `A == B`
end
