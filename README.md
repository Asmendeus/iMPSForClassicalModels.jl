# iMPSForClassicalModels.jl

`iMPSForClassicalModel.jl` is a package based on `TensorKit.jl` for solving 2D classical statistical models under thermodynamic limit.

!WARNING: The package is still under rapid development, whose current version `v0.1` is a demo version. Some features on "top" level are available, but the most efficient implementation is not provided. These problems are expected to be improved in version `v0.2`

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
