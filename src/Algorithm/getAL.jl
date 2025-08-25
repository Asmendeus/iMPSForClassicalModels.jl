"""
    getAL(AC::LocalTensor{R}, C::BondTensor) -> AL::LocalTensor{R}

Generate left-canonical tensor `AL` by center tensor `AC` and center bond tensor `C`, where `AC = AL C`.
"""
function getAL(AC::LocalTensor{R}, C::BondTensor) where R
    R ≥ 3 || throw(ArgumentError("Setting `BondTensor` as a center tensor is not supported"))
    UL_AC, _ = leftorth(AC; alg=Polar())
    UL_C, _ = leftorth(C; alg=Polar())
    return UL_AC * UL_C'
end
