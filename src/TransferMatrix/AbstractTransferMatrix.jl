"""
     abstract type AbstractTransferMatrix{L}

Wrapper type for transfer matrix.
"""
abstract type AbstractTransferMatrix{L} end

const AbstractMPSTransferMatrix = AbstractTransferMatrix{3}
const AbstractMPOTransferMatrix = AbstractTransferMatrix{4}
const AbstractMPSOrMPOTransferMatrix = Union{AbstractMPSTransferMatrix, AbstractMPOTransferMatrix}