# SparseMatricesForParallelAssembly
In this repository I investigate different sparse matrix formats which can be used for the assembly process in FVM (or FEM). Using Julia.

We divide the "assembly" process into 3 steps:

## 1. The true assembly
Which format performs best while building up the matrix?

## 2. Conversion to CSC
Which format / which algorithm performs best for converting the sum of multiple matrices to a CSC matrix?

## 3. Assembly in CSC
For example when using the Newton method to solve a nonlinear PDE, we only use 1. to find the structure of the matrix, but we directly change the entries of the CSC matrix
