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
        alg_grad::GradientAlgorithm
        alg_grad_mutiply::GradientAlgorithm
        alg_eig_mutiply::EigenAlgorithm
    end

ViTEBD algorithm struct storing parameters.
"""
struct ViTEBD <: FixedBondaryAlgorithm
    alg_grad::GradientAlgorithm
    alg_grad_mutiply::GradientAlgorithm
    alg_eig_mutiply::EigenAlgorithm

    function ViTEBD(alg_grad::GradientAlgorithm=SimpleIterator(;tol=Defaults.tol_low),
                    alg_grad_mutiply::GradientAlgorithm=Defaults.alg_grad,
                    alg_eig_mutiply::EigenAlgorithm=Defaults.alg_eig)
        return new(alg_grad, alg_grad_mutiply, alg_eig_mutiply)
    end
end
