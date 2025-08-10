function isSparse(ρ::SparseMPO; dense::Float64=Defaults.dense)
end

function isHermitian(ρ::SparseMPO; tol::Float64=Defaults.tol_norm)
    # conjugate + switch up and down
end