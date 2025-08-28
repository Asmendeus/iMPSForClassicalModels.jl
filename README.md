# iMPSForClassicalModels.jl

`iMPSForClassicalModel.jl` is a package based on `TensorKit.jl` for solving 2D classical statistical models under thermodynamic limit.

## Problems can be solved

`iMPSForClassicalModel.jl` provides a range of features for MPS/MPO with canonical form and implementations of some algorithms (including `ViTEBD`and `VUMPS`). Here's what it can solve

* Solve 2D classical model under thermodynamic limit (including thermodynamic quantities, 1-site physical quantity expectations, and 2-site correlation functions);
* Fastly implement any functions based on iMPS/iMPO;
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

## About Development and Maintenance

An iMPS package for quantum systems is planned to be developed.

A development plan is to incorporate the functionality implemented by `iMPSForClassicalModels.jl` into it, at which point maintenance of this package will cease.
