# SparseMatricesForParallelAssembly
In this repository I investigate different sparse matrix formats which can be used for the assembly process in FVM (or FEM). Using Julia.

## Problem statement

We divide the "assembly" process into 3 steps:

### 1. The true assembly
Which format performs best while building up the matrix?

### 2. Conversion to CSC
Which format / which algorithm performs best for converting the sum of multiple matrices to a CSC matrix?

### 3. Assembly in CSC
For example when using the Newton method to solve a nonlinear PDE, we only use 1. to find the structure of the matrix, but we directly change the entries of the CSC matrix

## General Idea

### 1.

### 2.

### 3.


## How to benchmark the code

Use `src/auto.jl`. 
The `path` is set automatically using `pwd()`.
The `pre_name` has to be set manually to a directory (which has to exist when running the code) where the data should be stored.
Then some 2d and 3d grids are created and the assembly, the conversion to CSC and the CSC-assembly are benchmarked. Using the `ns` arrays you can change which grids you want to test.

**Caution: The preparation of the grids (i.e. conversion to graph and partitioning) are very slow for fine/large grids, especially in 3d.**











