module Defaults
    datatype = ComplexF64

    dense = 0.1         # sparse matrix is matrix with dense degree smaller than it

    tol = 1e-14         # iteration tolerance
    maxiter = 200       # iteration step limit

    tol_norm = 1e-14    # norm(A - B): tolerance of equation `A == B`
end
