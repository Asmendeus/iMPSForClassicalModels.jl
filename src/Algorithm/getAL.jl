"""
    getAL(AC::LocalTensor{R}, C::BondTensor) -> AL::LocalTensor{R}

Generate left-canonical tensor `AL` by center tensor `AC` and center bond tensor `C`, where `AC = AL C`.
"""
function getAL(AC::LocalTensor{R}, C::BondTensor) where R
    R ≥ 3 || throw(ArgumentError("Setting `BondTensor` as a center tensor is not supported"))
    UL_AC, _ = leftorth(AC, Tuple(1:R-1), (R,); alg=Polar())
    UL_C, _ = leftorth(C, (1,), (2,); alg=Polar())
    return UL_AC * UL_C'
end
