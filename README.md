# iMPSForClassicalModels.jl

`iMPSForClassicalModel.jl` is a package based on `TensorKit.jl` for solving 2D classical statistical models under thermodynamic limit.

## Problems can be solved

`iMPSForClassicalModel.jl` provides a range of features for MPS/MPO with canonical form and implementations of some algorithms (including `ViTEBD`and `VUMPS`). Here's what it can solve

* Solve 2D classical model under thermodynamic limit (including thermodynamic quantities, 1-site physical quantity expectations, and 2-site correlation functions);
* Fastly implement any functions based on MPS/MPO with canonical form;
* Solve the ground state of a 1D quantum model under thermodynamic limit based on the principle of imaginary time evolution (controlling global quantum numbers is not supported).

## Package structure

```
### top ###
MPS/		# iMPS as a variational boundary
MPO/		# iMPO as a variational boundary
SparseMPO/	# iMPO corresponding to partition function
LocalImpurity/	# impurity tensor corresponding to physical quantity
Algorithm/	# algorithms based on MPS/MPO

### middle ###
Environment/	# local structure consisting of a few tensors, corresponding to fixed point equations
Method/		# methods targeting basic tensor or local structure

### ground ###
TensorWrapper/	# basic tensor type: classify `TensorMap`
```

## About tensor for classical model

See [Review.md](docs/Review.md)

## About infinite MPS

See [Review.md](docs/Review.md)

## Tutorial

See [Tutorial.md](docs/Tutorial.md)

## About next version

In `v0.2.0`, the following features are expected to be added

- Variational optimization implementation based on gradient algorithm
- Allowing `leftFixedPoint` and `rightFixedPoint` to specify solving which fixed point equation
- Optimization of partial code logic
- Examples of classical Clock and XY model
- A more accurate and comprehensive tutorial

## About the reference package - `FiniteMPS.jl`

`iMPSForClassicalModel.jl` refers to `FiniteMPS.jl v1.6.1` during the framework design process. If you are a user of it, note the following differences between data types and functions of the same or similar name.

For `TensorWrapper`,

- `MPSTensor{R}` in `FiniteMPS.jl` is renamed  `LocalTensor{R}`, and `const MPSTensor = LocalTensor{R}`;
- `leftorth` and `rightorth` return instances of type `AbstractTensorWrapper` but not `AbstractTensorMap`;

For `Method`,

- `pushleft` and `pushright` are pure functions but not mutating functions, used to find fixed point of various environments in infinite system;

For `iMPS`, similar to `MPS` in `FiniteMPS.jl`

- `canonicalize!` support only `si` as the orthogonal center, not `siL` and `siR`;

For `iMPO`, similar to `MPO` in `FiniteMPS.jl`

- `identityInfiniteMPO` generates a direct product state as an iMPO, made up of onsite identity matrices;

For `SparseMPO`,

**Note**: `SparseMPO` is an iMPO type used to store local tensors of the classical system's partition function. Although it is written as "SparseMPO", it is not always sparse for a classical system, which is just a naming convention inherited from `FiniteMPS.jl`. For a classical model with discrete degrees of freedom (like Ising, Clock), the local tensor corresponding to the partition function is usually dense, with dense degree equal to about 1. However, for a XY-like model with continuous degrees of freedom, the local tensors are usually sparse after discretization. For example, dense degree of classical XY model with truncation dimension 5 (physical dimension d = 11) on square lattice approximates to 0.06, which is typically sparse.

- `SparseMPO`'s field `A` is matrix of `W` rows and `L` columns, but not a vector of length `L`;
- Eltype of `SparseMPO.A` is `MPOTensor (alias for LocalTensor{4} where T)`, but not a sparse tensor type.
