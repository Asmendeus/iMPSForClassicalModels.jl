module Defaults
    datatype = ComplexF64

    tol = 1e-10          # iteration tolerance
    maxiter = 400        # iteration step limit

    tol_norm = 1e-12     # norm(A - B): tolerance of equation `A == B`
end
