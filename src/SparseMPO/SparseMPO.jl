"""
    struct SparseMPO{W, L} <: AbstractInfiniteMPO{L}
        O::Matrix{MPOTensor}
    end

!Warning: "SparseMPO" is a MPO type used to local tensors of the classical system's partition function. Although it is written as "SparseMPO", it is not always sparse for classical systems.
For a classical model with discrete degrees of freedom (like Ising, Clock), the local tensor corresponding to the partition function is usually dense, with dense degree equal to about 1.
However, for a XY-like model with continuous degrees of freedom, the local tensors are usually sparse after discretization.
For example, dense degree of classical XY model with truncation dimension 5 (physical dimension d = 11) on square lattice approximates to 0.06, which is typically sparse.

# TODO
Optimization for solving the ground state of dense mpo
"""
