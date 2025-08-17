"""
    getAR(AC::LocalTensor{R}, C::BondTensor) -> AR::LocalTensor{R}

Generate right-canonical tensor `AR` by center tensor `AC` and center bond tensor `C`, where `AC = C AR`.
"""
function getAR(AC::LocalTensor{R}, C::BondTensor) where R
    R ≥ 3 || throw(ArgumentError("Setting `BondTensor` as a center tensor is not supported"))
    _, UR_AC = rightorth(AC, (1,), Tuple(2:R); alg=Polar())
    _, UR_C =rightorth(C, (1,), (2,); alg=Polar())
    return BondTensor(UR_C.A') * UR_AC
end
