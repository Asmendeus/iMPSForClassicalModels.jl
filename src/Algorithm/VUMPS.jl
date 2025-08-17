"""
    struct VUMPS
        alg_eig::EigenAlgorithm
        tol_eig::Float64

        alg_grad::GradientAlgorithm
        tol_grad::Float64
    end

VUMPS algorithm struct storing parameters.
"""
struct VUMPS
    alg_eig::EigenAlgorithm
    tol_eig::Float64

    alg_grad::GradientAlgorithm
    tol_grad::Float64
end
