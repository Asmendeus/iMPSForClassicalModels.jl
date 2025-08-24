abstract type FixedBondaryAlgorithm end

"""
    struct VUMPS
        alg_eig::EigenAlgorithm
        tol_eig::Float64

        alg_grad::GradientAlgorithm
        tol_grad::Float64
    end

VUMPS algorithm struct storing parameters.
"""
struct VUMPS <: FixedBondaryAlgorithm
    alg_eig::EigenAlgorithm
    tol_eig::Float64

    alg_grad::GradientAlgorithm
    tol_grad::Float64
end

"""
    struct ViTEBD <: FixedBondaryAlgorithm
        alg_eig::EigenAlgorithm
        tol_eig::Float64

        alg_grad::GradientAlgorithm
        tol_grad::Float64

        maxlayer::Int64
    end

ViTEBD algorithm struct storing parameters.
"""
struct ViTEBD <: FixedBondaryAlgorithm
    alg_eig::EigenAlgorithm
    tol_eig::Float64

    alg_grad::GradientAlgorithm
    tol_grad::Float64

    width_max::Int64
end
