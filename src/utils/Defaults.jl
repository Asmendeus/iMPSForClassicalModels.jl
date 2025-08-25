module Defaults
    using iMPSForClassicalModels
    using KrylovKit

    datatype = ComplexF64

    dense = 0.1         # sparse matrix is matrix with dense degree smaller than it

    tol_high = 1e-14        # high-accuracy iteration tolerance
    tol = 1e-12             # normal-accuracy iteration tolerance
    tol_low = 1e-10         # low-accuracy iteration tolerance

    maxiter = 100       # iteration step limit

    alg_eig = KrylovKit.Arnoldi()           # algorithm for solving maximum/minimum eigenvalue problem
    alg_grad = iMPSForClassicalModels.SimpleIterator(;tol=tol, maxiter=maxiter)  # algorithm for solving gradient problem
end
