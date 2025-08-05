abstract type AbstractMPSTensor end

"""
       |
    —— A ——
"""
struct MPSTensor <: AbstractMPSTensor
end

"""
       |
    —— AL ——
"""
struct LeftIsometricTensor <: AbstractMPSTensor
end

"""
       |
    —— AR ——
"""
struct RightIsometricTensor <: AbstractMPSTensor
end