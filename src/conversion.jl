include(path*"preparatory.jl")

"""
`function CSC_RLNK_si_oc_ps_dz(As::Vector{SparseMatrixLNK{Float64, Int32}}, nr, s, nt, depth)`

Almost exactly `CSC_RLNK_si_oc_ps_dz`. Only difference, the zero entries of the LNKs are filtered out.

"""
function CSC_RLNK_si_oc_ps_dz(
	As::Vector{SparseMatrixLNK{Float64, Int32}}, 
	nr, s, nt, depth
	)

	nnz = sum([As[i].nnz for i=1:depth*nt+1]) #you could also subtract the diagonal entries from shared columns, since those are definitely double
	indptr = zeros(Int32, As[1].m+1)
	indices = zeros(Int32, nnz) #sum(As.nnz))
	data = zeros(Float64, nnz) #sum(As.nnz))
	ctr = 1
	eqctr = 0
	
	for j=1:As[1].m
		indptr[j] = ctr
		regionctr = 1
		jc = 0
		nrr = @view nr[:,j] 
		for region in @view nr[:,j] #nrr #[:,j]
			if region > 0
				k = s[Int(ceil(region/nt)), j]
				if regionctr == 1
					while k>0
						if As[region].nzval[k] != 0.0
							indices[ctr] = As[region].rowval[k]
							data[ctr]    = As[region].nzval[k]
							
							for jcc=1:jc
								if indices[ctr-jcc] > indices[ctr-jcc+1]
									tmp_i = indices[ctr-jcc+1]
									tmp_d = data[ctr-jcc+1]
									indices[ctr-jcc+1] = indices[ctr-jcc]
									data[ctr-jcc+1]    = data[ctr-jcc]
									
									indices[ctr-jcc] = tmp_i
									data[ctr-jcc]    = tmp_d
								else
									break
								end
							end
							
							ctr += 1
							jc += 1
						end
						k = As[region].colptr[k]
					end
				else
					while k>0
						if As[region].nzval[k] != 0.0
							indices[ctr] = As[region].rowval[k]
							data[ctr]    = As[region].nzval[k]
							
							for jcc=1:jc
								if indices[ctr-jcc] > indices[ctr-jcc+1]
									tmp_i = indices[ctr-jcc+1]
									tmp_d = data[ctr-jcc+1]
									indices[ctr-jcc+1] = indices[ctr-jcc]
									data[ctr-jcc+1]    = data[ctr-jcc]
									
									indices[ctr-jcc] = tmp_i
									data[ctr-jcc]    = tmp_d
								elseif indices[ctr-jcc] == indices[ctr-jcc+1]
									data[ctr-jcc] += data[ctr-jcc+1]
									eqctr += 1
									
									for jccc=1:jcc
										indices[ctr-jcc+jccc] = indices[ctr-jcc+jccc+1]
										data[ctr-jcc+jccc]    = data[ctr-jcc+jccc+1]
									end
									
									ctr -= 1
									jc  -= 1
									
									break
								else
									break
								end
							end
							
							ctr += 1
							jc += 1
						end
						k = As[region].colptr[k]
					end
					
				end
				regionctr += 1
			end
		end
		
	end

	indptr[end] = ctr

	#resize!(indices, ctr-1)
	#resize!(data, ctr-1)
	#SparseArrays.SparseMatrixCSC(As[1].m, As[1].m, indptr, indices, data)
	SparseArrays.SparseMatrixCSC(
		As[1].m, As[1].m, indptr, indices[1:ctr-1], data[1:ctr-1]
	)
	
end


"""
`function CSC_RLNK_si_oc_ps(As::Vector{SparseMatrixLNK{Float64, Int32}}, nr, s, nt, depth)`

This function creates a CSC matrix filled with the sum of the LNK matrices that are entered as elements of `As`. BUT: It's very different from `CSC_LNKs_s!`, because the submatrices are reduced in size.
We know beforehand into which columns of the 'global' matrix each submatrix can write, thus the submatrices only have a reduced number of columns. `s` which is `sortednodesperthread` from `get_nnnts_and_sortednodesperthread_and_noderegs_from_cellregs_ps` in `preparatory.jl`.
`nr` is the noderegions, `nt` the number of threads, `nt` is the number of threads and `depth` is the number of separator-partitonings (successively).

"""
function CSC_RLNK_si_oc_ps(
	As::Vector{SparseMatrixLNK{Float64, Int32}}, 
	nr, s, nt, depth
	)

	nnz = sum([As[i].nnz for i=1:depth*nt+1]) #you could also subtract the diagonal entries from shared columns, since those are definitely double
	indptr = zeros(Int32, As[1].m+1)
	indices = zeros(Int32, nnz) #sum(As.nnz))
	data = zeros(Float64, nnz) #sum(As.nnz))
	ctr = 1
	eqctr = 0
	
	for j=1:As[1].m
		indptr[j] = ctr
		regionctr = 1
		jc = 0
		nrr = @view nr[:,j] 
		for region in @view nr[:,j] #nrr #[:,j]
			if region > 0
				k = s[Int(ceil(region/nt)), j]
				if regionctr == 1
					while k>0
						indices[ctr] = As[region].rowval[k]
						data[ctr]    = As[region].nzval[k]
						
						for jcc=1:jc
							if indices[ctr-jcc] > indices[ctr-jcc+1]
								tmp_i = indices[ctr-jcc+1]
								tmp_d = data[ctr-jcc+1]
								indices[ctr-jcc+1] = indices[ctr-jcc]
								data[ctr-jcc+1]    = data[ctr-jcc]
								
								indices[ctr-jcc] = tmp_i
								data[ctr-jcc]    = tmp_d
							else
								break
							end
						end
						
						ctr += 1
						jc += 1
						k = As[region].colptr[k]
					end
				else
					while k>0
						indices[ctr] = As[region].rowval[k]
						data[ctr]    = As[region].nzval[k]
						
						for jcc=1:jc
							if indices[ctr-jcc] > indices[ctr-jcc+1]
								tmp_i = indices[ctr-jcc+1]
								tmp_d = data[ctr-jcc+1]
								indices[ctr-jcc+1] = indices[ctr-jcc]
								data[ctr-jcc+1]    = data[ctr-jcc]
								
								indices[ctr-jcc] = tmp_i
								data[ctr-jcc]    = tmp_d
							elseif indices[ctr-jcc] == indices[ctr-jcc+1]
								data[ctr-jcc] += data[ctr-jcc+1]
								eqctr += 1
								
								for jccc=1:jcc
									indices[ctr-jcc+jccc] = indices[ctr-jcc+jccc+1]
									data[ctr-jcc+jccc]    = data[ctr-jcc+jccc+1]
								end
								
								ctr -= 1
								jc  -= 1
								
								break
							else
								break
							end
						end
						
						ctr += 1
						jc += 1
						k = As[region].colptr[k]
					end
					
				end
				regionctr += 1
			end
		end
		
	end

	indptr[end] = ctr
	
	#resize!(indices, ctr-1)
	#resize!(data, ctr-1)
	#SparseArrays.SparseMatrixCSC(As[1].m, As[1].m, indptr, indices, data)
	SparseArrays.SparseMatrixCSC(
		As[1].m, As[1].m, indptr, indices[1:ctr-1], data[1:ctr-1]
	)
	
end




"""
`function CSC_LNK_si(lnk::SparseMatrixLNK{Tv, Ti}) where {Tv, Ti <: Integer}`

Turns the LNK matrix `lnk` into a CSC matrix.
Stolen from `https://github.com/j-fu/ExtendableSparse.jl/blob/master/src/matrix/sparsematrixlnk.jl`, `https://github.com/j-fu/ExtendableSparse.jl/blob/master/src/matrix/sparsematrixlnk.jl`.. 
Only modification is, there is no CSC matrix on which the LNK matrix added, this is a speed up.
"""
function CSC_LNK_si(lnk::SparseMatrixLNK{Tv, Ti}) where {Tv, Ti <: Integer}
    csc = spzeros(lnk.m, lnk.n)
    # overallocate arrays in order to avoid
    # presumably slower push!
    xnnz = nnz(lnk)
    colptr = Vector{Ti}(undef, csc.n + 1)
    rowval = Vector{Ti}(undef, xnnz)
    nzval = Vector{Tv}(undef, xnnz)

	l_lnk_col = 0
        
    inz = 1 # counts the nonzero entries in the new matrix
	for j = 1:(csc.n)
		jc = 0
        # Copy extension entries into col and sort them
        k = j
        colptr[j] = inz
        while k > 0
            if lnk.rowval[k] > 0
                l_lnk_col += 1
				rowval[l_lnk_col] = lnk.rowval[k]
				nzval[l_lnk_col]  = lnk.nzval[k]
				
                #col[l_lnk_col] = ColEntry(lnk.rowval[k], lnk.nzval[k])

				for jcc=1:jc
					if rowval[l_lnk_col-jcc] > rowval[l_lnk_col-jcc+1]
						tmp_r = rowval[l_lnk_col-jcc+1]
						tmp_v = nzval[l_lnk_col-jcc+1]
						rowval[l_lnk_col-jcc+1] = rowval[l_lnk_col-jcc]
						nzval[l_lnk_col-jcc+1]  = nzval[l_lnk_col-jcc]
						rowval[l_lnk_col-jcc] = tmp_r
						nzval[l_lnk_col-jcc]  = tmp_v
					else
						break
					end
						
				end
				inz += 1
				jc += 1
            end
            k = lnk.colptr[k]
        end
    end
    colptr[csc.n + 1] = inz
    # Julia 1.7 wants this correct
    resize!(rowval, inz - 1)
    resize!(nzval, inz - 1)
    SparseMatrixCSC(csc.m, csc.n, colptr, rowval, nzval)
end

"""
`function add_LNKs!(A,B)`

A += B for LNK matrices.
"""
function add_LNKs!(A,B) #A += B
	n = B.n #number of columns #size(B)[2]
	val = B.nzval
	col = B.colptr
	row = B.rowval

	j = 1
	id = 1
	ctr = 0
	while j <= n
		if row[id] > 0
			A[row[id], j] += val[id]
		end

		if col[id] == 0
			j += 1
			id = j
		else
			id = col[id]
		end
		if j > n
			break
		end
	end
end

"""
`function CSC_LNKs_s!(As::Vector{SparseMatrixLNK{Tv, Ti}}) where {Tv, Ti <: Integer}`

Creates a CSC matrix that contains the sum of the LNK matrices that are entered as elements of `As`.
In other words:
`As = [A1, A2, A3, A4]; A1 += A2; A1 += A2; A1 += A3; A1 += A4; return CSC(A1)`

"""
function CSC_LNKs_s!(As::Vector{SparseMatrixLNK{Tv, Ti}}) where {Tv, Ti <: Integer}
	A = CSC_LNK_si(As[1])
	for i=2:length(As)
		A += CSC_LNK_si(As[i])
	end
	A
end


### Landfill:

function nnzincolumn(A, j)
	k = j
	ctr = 1
	while k > 0
		k = A.colptr[k]
		ctr += 1
	end
	ctr-1
end

function get_nnzc(As, nt, nr, s)
	nnzc = zeros(Int32, As[1].m)
	
	for j=1:As[1].m
		for region in @view nr[:,j] #nrr #[:,j]
			if region > 0
				k = s[Int(ceil(region/nt)), j]
				while k>0
					k = As[region].colptr[k]
					nnzc[j] += 1
				end
			end
		end
	end
	
	nnzc
end
