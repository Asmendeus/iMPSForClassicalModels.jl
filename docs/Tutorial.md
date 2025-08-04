# Tutorial

## 0. Pre-knowledge

### Partition function

For any two-dimensional discrete classical statistical model under thermodynamic limit with a Hamiltonian $H$, the partition function $Z = \mathrm{e}^{-\beta H}$ can always be expressed as an infinite tensor network structure (perhaps after some approximation or transformation) as follow

![PartitionFunction](fig/PartitionFunction.png "Partition Function")

The value of the partition function $Z$ is the result of the contraction of the tensor network.

To solve it, we can consider a variational boundary iMPS at its lower  boundary at infinity as follow

![BoundaryMPS](fig/BoundaryMPS.png "Boundary MPS")

Next, we push the boundary MPS up layer by layer to contract. If we limit its bond dimension to $D$, the boundary iMPS will eventually converge to a fixed point. The fixed point is the solution of the following eigen equation,

![EigenEquation](fig/EigenEquation.png "Eigen Equation")

where $\lambda$ is the dominant eigenvalue. Such a maximum eigenvalue problem can be solved by `VUMPS` or  `iTEBD` algorithms.

In addition, given the translational symmetry of iMPS, we are actually solving the following equation in the iMPS environment

![EigenEquation2](fig/EigenEquation2.png "Local Eigen Equation")

where $\lambda_i$ is the dominant eigenvalue, which is local partition function as well. The total partition function is the product of all local partition functions, i.e.,

$$
Z=\prod_{i} \lambda_i
$$

### 1-site physical quantity

### 2-site correlation function

### Two generalizations

# 1.
