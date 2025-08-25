"""
    getAR(AC::LocalTensor{R}, C::BondTensor) -> AR::LocalTensor{R}

Generate right-canonical tensor `AR` by center tensor `AC` and center bond tensor `C`, where `AC = C AR`.
"""
function getAR(AC::LocalTensor{R}, C::BondTensor) where R
    R ≥ 3 || throw(ArgumentError("Setting `BondTensor` as a center tensor is not supported"))
    _, UR_AC = rightorth(AC; alg=Polar())
    _, UR_C =rightorth(C; alg=Polar())
    return UR_C' * UR_AC
end
