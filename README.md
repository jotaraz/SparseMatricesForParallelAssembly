# SparseMatricesForParallelAssembly
In this repository I investigate different sparse matrix formats which can be used for the assembly process in FVM (or FEM). Using Julia.

For now I have looked at LNK (SparseMatrixLNK from https://github.com/j-fu/ExtendableSparse.jl/tree/master) and a dictionary based format (SparseMatrixDict from https://github.com/masuday/SparseMatrixDicts.jl/tree/master).
For those two formats I have looked at twoi different parts of the "assembly" process.
## 1. The true assembly
Which format performs best while building up the matrix?

## 2. Conversion to CSC
Which format / which algorithm performs best for converting the sum of multiple matrices to a CSC matrix?

