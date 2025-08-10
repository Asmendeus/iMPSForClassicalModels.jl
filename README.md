# iMPSForClassicalModels.jl

`iMPSForClassicalModel.jl` is a package based on `TensorKit.jl` for solving 2D classical statistical models under thermodynamic limit.

## Problems can be solved

`iMPSForClassicalModel.jl` provides a range of features for MPS/MPO with canonical form and implementations of some algorithms (including `ViTEBD`and `VUMPS`). Here's what it can solve

* Solve 2D classical model under thermodynamic limit (including thermodynamic quantities, 1-site physical quantity expectations, and 2-site correlation functions);
* Fastly implement any functions based on MPS/MPO with canonical form;
* Solve the ground state of a 1D quantum model under thermodynamic limit based on the principle of imaginary time evolution (controlling global quantum numbers is not supported).

## About tensor for classical model

See [Review.md](docs/Review.md)

## About infinite MPS

See [Review.md](docs/Review.md)

## Tutorial

See [Tutorial.md](docs/Tutorial.md)
