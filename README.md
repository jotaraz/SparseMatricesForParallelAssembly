# SparseMatricesForParallelAssembly
In this repository I investigate different sparse matrix formats which can be used for the assembly process in FVM (or FEM). Using Julia.

## Table of contents:
- Problem statement
- General idea
- How to benchmark the code
- To-Do

## Problem statement

We divide the "assembly" process into 3 steps:

### 1. The true assembly
Which format performs best while building up the matrix?

### 2. Conversion to CSC
Which format / which algorithm performs best for converting the sum of multiple matrices to a CSC matrix?

### 3. Assembly in CSC
For example when using the Newton method to solve a nonlinear PDE, we only use 1. to find the structure of the matrix, but we directly change the entries of the CSC matrix

There is a 0-th step:

### 0. Preparation
The grid needs to be partitioned...

## General idea

### 1. The true assembly
Our approach is based on linked-lists-matrices (LNK) using the [ExtendableSparse.jl](https://github.com/j-fu/ExtendableSparse.jl) implementation. 
The general idea is the following: Consider 4 threads, then the cells in your grid are partitioned into 4 regions with a separator, such that no node is in two partitions (but it could be in a partition and in the separator). After partitioning the cells, for each node we know in which region it is (and thus, which thread will access it). Then, we create 4+1 LNK (reduced) matrices, such that each matrix has as many columns as there are nodes in the respective region (or in the separator). We also have a mapping from indices in the reduced matrix to indices in the (big) matrix.
Then, each thread iterates over the cells in the respective region and writes into the reduced respective matrix. Afterwards, we iterate over the separator.
In the end we have 5 LNK reduced matrices which are all smaller then the original matrix.
This format is called RLNKs.

### 2. Conversion to CSC
In the last paragraph we explained we needed a mapping from indices in the submatrix to indices in the (big) matrix. We also need the reverse mapping, such that we iterate over the global node indices, get the column in the submatrix, put it in the [CSC](https://docs.julialang.org/en/v1/stdlib/SparseArrays/) `rowval`- and `nzval`-arrays. For the nodes that are also accessed by the separator, we need to access the separator matrix as well.

### 3. Assembly in CSC
**This is the current slow algorithm.**
After the initial assembly and the conversion we have the CSC matrix and we do not longer need the RLNK matrices. But, for solving nonlinear PDEs, the Newton method is used where the values in the CSC matrix are changed. Most of the times, only the values are changed, not the structure.
But somtimes there is a value that could be $\neq 0$ (structurally) but happens to be 0 in the initial assembly. Then this entry is not inlcuded in the CSC matrix. But in the next Newton step, this value could actually be $\neq 0$. For such cases, we have so-called backup reduced LNK matrices. After iterating over all cells, the CSC matrix has changed entries and we probably have some entries in the backup matrix, those need to be summed together.


## How to benchmark the code

Use [src/auto.jl](https://github.com/jotaraz/SparseMatricesForParallelAssembly/blob/main/src/auto.jl). 
The `path` is set automatically using `pwd()`.
The `pre_name` has to be set manually to a directory (which has to exist when running the code) where the data should be stored.
Then some 2d and 3d grids are created and the assembly, the conversion to CSC and the CSC-assembly are benchmarked. Using the `ns` arrays you can change which grids you want to test.

Using the [Pluto](https://github.com/fonsp/Pluto.jl) notebook [src/plot_benchmark_data.jl](https://github.com/jotaraz/SparseMatricesForParallelAssembly/blob/main/src/plot_benchmark_data.jl) the benchmark data I already collected can be visulaized.

**Caution: The preparation of the grids (i.e. conversion to graph and partitioning) are very slow for fine/large grids, especially in 3d.**

## To-Do

- Optimize CSC-assembly
- Optimize preparation
- ILUZero
- Boundary treatment
- Implementation (in VoronoiFVM.jl)









