using SparseArrays
using ExtendableSparse
using Base.Threads
using ThreadPinning
pinthreads(:cores)

nt = nthreads()


function array_ci(x)
	if length(x)==0
		return "length 0 "
	else
		return "length $(length(x)), max $(maximum(abs.(x))) "
	end
end

"""
`function pe_s(lnk::SparseMatrixLNK{Tv, Ti}, csc::SparseMatrixCSC)::SparseMatrixCSC where {Tv, Ti <: Integer}`

This function is taken from Jürgen Fuhrmann's `https://github.com/j-fu/ExtendableSparse.jl/blob/master/src/matrix/sparsematrixlnk.jl` the function there is called `function Base.:+(lnk::SparseMatrixLNK{Tv, Ti}, csc::SparseMatrixCSC)::SparseMatrixCSC where {Tv, Ti <: Integer}`
"""
function pe_s(lnk::SparseMatrixLNK{Tv, Ti},
                 csc::SparseMatrixCSC)::SparseMatrixCSC where {Tv, Ti <: Integer}
    @assert(csc.m==lnk.m)
    @assert(csc.n==lnk.n)

    # overallocate arrays in order to avoid
    # presumably slower push!
    xnnz = nnz(csc) + nnz(lnk)
    colptr = Vector{Ti}(undef, csc.n + 1)
    rowval = Vector{Ti}(undef, xnnz)
    nzval = Vector{Tv}(undef, xnnz)

    # Detect the maximum column length of lnk
    lnk_maxcol = 0
    for j = 1:(csc.n)
        lcol = zero(Ti)
        k = j
        while k > 0
            lcol += 1
            k = lnk.colptr[k]
        end
        lnk_maxcol = max(lcol, lnk_maxcol)
    end

    # pre-allocate column  data
    col = [ColEntry{Tv, Ti}(0, zero(Tv)) for i = 1:lnk_maxcol]

    inz = 1 # counts the nonzero entries in the new matrix

    in_csc_col(jcsc, j) = (nnz(csc) > zero(Ti)) && (jcsc < csc.colptr[j + 1])

    in_lnk_col(jlnk, l_lnk_col) = (jlnk <= l_lnk_col)

    # loop over all columns
    for j = 1:(csc.n)
        # Copy extension entries into col and sort them
        k = j
        l_lnk_col = 0
        while k > 0
            if lnk.rowval[k] > 0
                l_lnk_col += 1
                col[l_lnk_col] = ColEntry(lnk.rowval[k], lnk.nzval[k])
            end
            k = lnk.colptr[k]
        end
        sort!(col, 1, l_lnk_col, Base.QuickSort, Base.Forward)

        # jointly sort lnk and csc entries  into new matrix data
        # this could be replaced in a more transparent manner by joint sorting:
        # make a joint array for csc and lnk col, sort them.
        # Will this be faster? 

        colptr[j] = inz
        jlnk = one(Ti) # counts the entries in col
        jcsc = csc.colptr[j]  # counts entries in csc

        while true
            if in_csc_col(jcsc, j) &&
               (in_lnk_col(jlnk, l_lnk_col) && csc.rowval[jcsc] < col[jlnk].rowval ||
                !in_lnk_col(jlnk, l_lnk_col))
                # Insert entries from csc into new structure
                rowval[inz] = csc.rowval[jcsc]
                nzval[inz] = csc.nzval[jcsc]
                jcsc += 1
                inz += 1
            elseif in_csc_col(jcsc, j) &&
                   (in_lnk_col(jlnk, l_lnk_col) && csc.rowval[jcsc] == col[jlnk].rowval)
                # Add up entries from csc and lnk
                rowval[inz] = csc.rowval[jcsc]
                nzval[inz] = csc.nzval[jcsc] + col[jlnk].nzval
                jcsc += 1
                inz += 1
                jlnk += 1
            elseif in_lnk_col(jlnk, l_lnk_col)
                # Insert entries from lnk res. col into new structure
                rowval[inz] = col[jlnk].rowval
                nzval[inz] = col[jlnk].nzval
                jlnk += 1
                inz += 1
            else
                break
            end
        end
    end
    colptr[csc.n + 1] = inz
    # Julia 1.7 wants this correct
    resize!(rowval, inz - 1)
    resize!(nzval, inz - 1)
    SparseMatrixCSC{Tv, Ti}(csc.m, csc.n, colptr, rowval, nzval)
end


struct ColEntry{Tv, Ti <: Integer}
    rowval::Ti
    nzval::Tv
end

# Comparison method for sorting
Base.isless(x::ColEntry, y::ColEntry) = (x.rowval < y.rowval)

function pe_p(lnk::SparseMatrixLNK{Tv, Ti},
                 csc::SparseMatrixCSC, nt)::SparseMatrixCSC where {Tv, Ti <: Integer}
    @assert(csc.m==lnk.m)
    @assert(csc.n==lnk.n)

	nc_nt = Int(csc.n/nt)

    
    # Detect the maximum column length of lnk
    lnk_maxcols = zeros(Ti, nt)
    xnnzs  = [csc.colptr[tid*nc_nt+1]-csc.colptr[(tid-1)*nc_nt+1] for tid=1:nt]
     
    @threads :static for tid=1:nt
    	inz = 1
    	for j=(tid-1)*nc_nt+1:tid*nc_nt
    		lcol = zero(Ti)
    		k = j
    		while k > 0
    			lcol += 1
            	k = lnk.colptr[k]
        	end
        	xnnzs[tid] += lcol
        	lnk_maxcols[tid] = max(lcol, lnk_maxcols[tid])
    	end
    end
	#=    
    for j = 1:(csc.n)
        lcol = zero(Ti)
        k = j
        while k > 0
            lcol += 1
            k = lnk.colptr[k]
        end
        lnk_maxcol = max(lcol, lnk_maxcol)
    end
	=#
	# overallocate arrays in order to avoid
    # presumably slower push!
    xnnzs[1]   = nnz(csc) + nnz(lnk)
    lengths    = nc_nt*ones(Ti, nt)
    lengths[1] = csc.n
    colptrs = [Vector{Ti}(undef, lengths[tid] + 1) for tid=1:nt]
    rowvals = [Vector{Ti}(undef, xnnzs[tid]) for tid=1:nt]
    nzvals  = [Vector{Tv}(undef, xnnzs[tid]) for tid=1:nt]

    # pre-allocate column  data
    cols = [[ColEntry{Tv, Ti}(0, zero(Tv)) for i = 1:lnk_maxcols[tid]] for tid=1:nt]

    #inzs = ones(Ti, nt) # counts the nonzero entries in the new matrix

    in_csc_col(jcsc, j) = (nnz(csc) > zero(Ti)) && (jcsc < csc.colptr[j + 1])

    in_lnk_col(jlnk, l_lnk_col) = (jlnk <= l_lnk_col)

    # loop over all columns
    @threads :static for tid=1:nt
    	inz = 1
    	for j=(tid-1)*nc_nt+1:tid*nc_nt
    		
    		# Copy extension entries into col and sort them
		    k = j
		    l_lnk_col = 0
		    while k > 0
		        if lnk.rowval[k] > 0
		            l_lnk_col += 1
		            cols[tid][l_lnk_col] = ColEntry(lnk.rowval[k], lnk.nzval[k])
		        end
		        k = lnk.colptr[k]
		    end
		    sort!(cols[tid], 1, l_lnk_col, Base.QuickSort, Base.Forward)

		    # jointly sort lnk and csc entries  into new matrix data
		    # this could be replaced in a more transparent manner by joint sorting:
		    # make a joint array for csc and lnk col, sort them.
		    # Will this be faster? 

		    colptrs[tid][j-(tid-1)*nc_nt] = inz
		    jlnk = one(Ti) # counts the entries in col
		    jcsc = csc.colptr[j]  # counts entries in csc

		    while true
		        if in_csc_col(jcsc, j) &&
		           (in_lnk_col(jlnk, l_lnk_col) && csc.rowval[jcsc] < cols[tid][jlnk].rowval ||
		            !in_lnk_col(jlnk, l_lnk_col))
		            # Insert entries from csc into new structure
		            rowvals[tid][inz] = csc.rowval[jcsc]
		            nzvals[tid][inz] = csc.nzval[jcsc]
		            jcsc += 1
		            inz += 1
		        elseif in_csc_col(jcsc, j) &&
		               (in_lnk_col(jlnk, l_lnk_col) && csc.rowval[jcsc] == cols[tid][jlnk].rowval)
		            # Add up entries from csc and lnk
		            rowvals[tid][inz] = csc.rowval[jcsc]
		            nzvals[tid][inz] = csc.nzval[jcsc] + cols[tid][jlnk].nzval
		            jcsc += 1
		            inz += 1
		            jlnk += 1
		        elseif in_lnk_col(jlnk, l_lnk_col)
		            # Insert entries from lnk res. col into new structure
		            rowvals[tid][inz] = cols[tid][jlnk].rowval
		            nzvals[tid][inz] = cols[tid][jlnk].nzval
		            jlnk += 1
		            inz += 1
		        else
		            break
		        end
		    end
    	end
    	colptrs[tid][nc_nt + 1] = inz
    end
    
    
    j0 = 0
    for tid=2:nt
    	j0 += colptrs[tid-1][nc_nt+1]-1
	    colptrs[1][(tid-1)*nc_nt+1:tid*nc_nt] = colptrs[tid][1:end-1].+j0
	    rowvals[1][j0+1:j0+colptrs[tid][end]-1] = view(rowvals[tid], 1:colptrs[tid][end]-1) #rowvals[tid][1:colptrs[tid][end]-1]
	    nzvals[1][j0+1:j0+colptrs[tid][end]-1]  = view(nzvals[tid], 1:colptrs[tid][end]-1) #nzvals[tid][1:colptrs[tid][end]-1]
    end
    colptrs[1][end] = j0+colptrs[nt][nc_nt+1]
    
    
    
    #colptr[csc.n + 1] = inz
    # Julia 1.7 wants this correct
    resize!(rowvals[1], colptrs[1][end] - 1)
    resize!(nzvals[1], colptrs[1][end] - 1)
    SparseMatrixCSC{Tv, Ti}(csc.m, csc.n, colptrs[1], rowvals[1], nzvals[1])
end
"""
`function compare_seq_to_parallel(num, n, k)`

1. Use the sequential `pe_s` and the parallel `pe_p` functions for `CSC + LNK` and compare the results.
2. Compare the time and allocations they take / make.

Random matrices are used. The number of threads used is the number Julia is started with (e.g. `julia -t 4`). `n` has to be an integer multiple of the number of threads. \\
num: Samplesize for the benchmark. \\
n: Both the LNK and the CSC matrix are square n x n matrices. \\
k: Number of changes in the matrix, k <= number of non-zero entries
"""
function compare_seq_to_parallel(num, n, k)
	times1 = zeros(num)
	times2 = zeros(num)
	A = SparseMatrixLNK{Float64, Int32}(n, n)
	B = SparseMatrixLNK{Float64, Int32}(n, n)

	for i=1:k
		A[rand(1:n), rand(1:n)] += 1.0
		B[rand(1:n), rand(1:n)] += 1.0
	end

	C = SparseMatrixCSC(A)

	C1 = pe_s(B, C)
	C2 = pe_p(B, C, nt)

	println("Length of (CSC_sequantial - CSC_parallel).nzval = $(length((C1-C2).nzval))")
	
	for i=1:num
		A = SparseMatrixLNK{Float64, Int32}(n, n)
		B = SparseMatrixLNK{Float64, Int32}(n, n)
		for i=1:k
			A[rand(1:n), rand(1:n)] += 1.0
			B[rand(1:n), rand(1:n)] += 1.0
		end
		C = SparseMatrixCSC(A)
		times1[i] = @elapsed pe_s(B, C)
		times2[i] = @elapsed pe_p(B, C, nt)
		GC.gc()
	end
	
	A = SparseMatrixLNK{Float64, Int32}(n, n)
	B = SparseMatrixLNK{Float64, Int32}(n, n)
	for i=1:k
		A[rand(1:n), rand(1:n)] += 1.0
		B[rand(1:n), rand(1:n)] += 1.0
	end
	C = SparseMatrixCSC(A)
	a1 = @allocated pe_s(B, C)
	a2 = @allocated pe_p(B, C, nt)
		
	println("Time  for sequential LNK+CSC: $(minimum(times1)) seconds")
	println("Time  for parallel   LNK+CSC: $(minimum(times2)) seconds")
	println("Allocs in sequential LNK+CSC: $a1 bytes")
	println("Allocs in parallel   LNK+CSC: $a2 bytes")
end
		
	
	
	
	
	


