struct ViTEBD
    alg_eig::Union{SimpleIteration, KrylovKit.KrylovAlgorithm}
    tol_eig::Float64

    alg_grad::Union{SimpleIteration}
    tol_grad::Float64

    maxlayer::Int64
end
