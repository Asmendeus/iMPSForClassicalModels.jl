"""
     abstract type AbstractTransferMatrix{L, R}

Wrapper type for transfer matrix.
"""
abstract type AbstractTransferMatrix{L, R} end

const AbstractMPSTransferMatrix{L} = AbstractTransferMatrix{L, 3}
const AbstractMPOTransferMatrix{L} = AbstractTransferMatrix{L, 4}
