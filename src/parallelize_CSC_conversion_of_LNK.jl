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


struct ColEntry{Tv, Ti <: Integer}
    rowval::Ti
    nzval::Tv
end

# Comparison method for sorting
Base.isless(x::ColEntry, y::ColEntry) = (x.rowval < y.rowval)



"""
y[1] = max(x[1,2,...,n]) \\
y[2] = max(x[n+1,n+2,...,2n]) \\
...
"""
function shortening(x, n)
	m = Int(length(x/n))
	y = zeros(typeof(x[1]), m)
	for i=1:m
		y[i] = maximum(x[(i-1)*n+1:i*n])
	end
	y	
end

function empirical_nt2()
	nt2 = Int(nt/2)
	if nt > 16
		nt2 = 8
	elseif nt <= 4
		nt2 = nt
	end
	nt2
end

# plusequal functions:
# - pe_s  : sequential
# - pe_p  : parallel
# - pe_p2 : parallel2


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




"""
`function pe_p(lnk::SparseMatrixLNK{Tv, Ti}, csc::SparseMatrixCSC, nt, nt2)::SparseMatrixCSC where {Tv, Ti <: Integer}`

First version of parallel plusequals.
We have multiple (shorter) loops over all threads. In comparison, `pe_p2` has only one and is slightly faster than `pe_p`.
The non-main loops in pe_p are very short, such that the multithreading overhead is large, thus only `nt2` threads are used there.
"""
function pe_p(lnk::SparseMatrixLNK{Tv, Ti},
                 csc::SparseMatrixCSC, nt, nt2)::SparseMatrixCSC where {Tv, Ti <: Integer}
    @assert(csc.m==lnk.m)
    @assert(csc.n==lnk.n)


	nc_nt = Int(csc.n/nt)
	
    
    # Detect the maximum column length of lnk
    lnk_maxcols = zeros(Ti, nt)
    xnnzs  = [csc.colptr[tid*nc_nt+1]-csc.colptr[(tid-1)*nc_nt+1] for tid=1:nt]
    
    fac = Int(nt/nt2)
	@threads :static for ltid=1:nt2
		for tid0=1:fac
			tid = (ltid-1)*fac+tid0
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
	end
    
	
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
    
    
    #times1[1+nt+1] = @elapsed begin
	j0 = 0
	for tid=2:nt
		j0 += colptrs[tid-1][nc_nt+1]-1
		colptrs[1][(tid-1)*nc_nt+1:tid*nc_nt] = colptrs[tid][1:end-1].+j0
	end
	colptrs[1][end] = j0+colptrs[nt][nc_nt+1]
	#end
	
	#times1[1+nt+2] = @elapsed begin
		#j0 = 0
	@threads :static for tid=2:nt
		j1 = colptrs[1][(tid-1)*nc_nt+1] #colptrs[tid-1][nc_nt+1]-1
		j2 = j1+colptrs[tid][end]-2 # colptrs[1][tid*nc_nt+1]-1
		rowvals[1][j1:j2] = view(rowvals[tid], 1:colptrs[tid][end]-1) #rowvals[tid][1:colptrs[tid][end]-1]
		nzvals[1][j1:j2]  = view(nzvals[tid], 1:colptrs[tid][end]-1) #nzvals[tid][1:colptrs[tid][end]-1]
	end
    
    
    
    
    #colptr[csc.n + 1] = inz
    # Julia 1.7 wants this correct
    resize!(rowvals[1], colptrs[1][end] - 1)
    resize!(nzvals[1], colptrs[1][end] - 1)
    SparseMatrixCSC{Tv, Ti}(csc.m, csc.n, colptrs[1], rowvals[1], nzvals[1])
end


"""
`function pe_p2(lnk::SparseMatrixLNK{Tv, Ti}, csc::SparseMatrixCSC, nt)::SparseMatrixCSC where {Tv, Ti <: Integer}`

2nd type of parallelized function of plusequals.
Here, we have only one big loop over all threads, this makes it a bit faster than `pe_p`.
"""
function pe_p2(lnk::SparseMatrixLNK{Tv, Ti},
                 csc::SparseMatrixCSC, nt)::SparseMatrixCSC where {Tv, Ti <: Integer}
    @assert(csc.m==lnk.m)
    @assert(csc.n==lnk.n)
    #@info csc.n
    
    stat = 0
    readys = zeros(Ti, nt)
    
    gxnnz = nnz(csc) + nnz(lnk)
    gcolptr = Vector{Ti}(undef, csc.n + 1)
    gcolptr[1] = 1
    growval = Vector{Ti}(undef, gxnnz)
    gnzval = Vector{Tv}(undef, gxnnz)
    nc_nt = Int(csc.n/nt)
    
    in_csc_col(jcsc, j) = (nnz(csc) > zero(Ti)) && (jcsc < csc.colptr[j + 1])
	in_lnk_col(jlnk, l_lnk_col) = (jlnk <= l_lnk_col)
	
	@threads :static for tid=1:nt
		xnnz  = csc.colptr[tid*nc_nt+1]-csc.colptr[(tid-1)*nc_nt+1]
			
		lnk_maxcol = zero(Ti)
			
		inz = 1
		for j=(tid-1)*nc_nt+1:tid*nc_nt
			lcol = zero(Ti)
			k = j
			while k > 0
				lcol += 1
			   	k = lnk.colptr[k]
		   	end
		   	xnnz += lcol
		   	lnk_maxcol = max(lcol, lnk_maxcol)
		end
		
		#colptr = Vector{Ti}(undef, nc_nt + 1) for tid=1:nt]
		rowval = Vector{Ti}(undef, xnnz) #for tid=1:nt]
		nzval  = Vector{Tv}(undef, xnnz) #for tid=1:nt]

		# pre-allocate column  data
		col = [ColEntry{Tv, Ti}(0, zero(Tv)) for i = 1:lnk_maxcol] # for tid=1:nt]

		#inzs = ones(Ti, nt) # counts the nonzero entries in the new matrix
		inz = 1
		
		for j=(tid-1)*nc_nt+1:tid*nc_nt
    		
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

			if j != (tid-1)*nc_nt+1
			    gcolptr[j] = inz
			end
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
    	
    	gcolptr[tid*nc_nt+1] = inz
		readys[tid] = 1
		#ready = sum(readys)
		#@info "while start $tid $stat $(sum(readys))" #tid, stat, ready
		while (stat != tid-1) || (sum(readys) != nt)
		end
		#@info "while done $tid $stat $(sum(readys))"
		
		j1 = gcolptr[(tid-1)*nc_nt + 1]
		j2 = j1+inz-2  #gcolptr[tid*nc_nt+1]-1
		
		
		if tid == 1
			#gcolptr[(tid-1)*nc_nt+2:end] .+= gcolptr[(tid-1)*nc_nt+1]-1
		else
			gcolptr[(tid-1)*nc_nt+2:end] .+= gcolptr[(tid-1)*nc_nt+1]-gcolptr[(tid-2)*nc_nt+1]
		end
		
		stat = tid
		
		growval[j1:j2] = view(rowval, 1:j2-j1+1) 
		gnzval[j1:j2]  = view(nzval, 1:j2-j1+1) 
		
    end
    
    resize!(growval, gcolptr[end] - 1)
    resize!(gnzval, gcolptr[end] - 1)
    SparseMatrixCSC{Tv, Ti}(csc.m, csc.n, gcolptr, growval, gnzval)
end


# timed plusequals functions
# these functions give a detailed view into which parts of the functions take how much time

"""
`function pe_s_time(lnk::SparseMatrixLNK{Tv, Ti}, csc::SparseMatrixCSC) where {Tv, Ti <: Integer}`

Timed version of `pe_s`.
"""
function pe_s_time(lnk::SparseMatrixLNK{Tv, Ti},
                 csc::SparseMatrixCSC) where {Tv, Ti <: Integer}
    @assert(csc.m==lnk.m)
    @assert(csc.n==lnk.n)
	times1 = zeros(3)
	times1[1] = @elapsed begin
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
	end
    # loop over all columns
    times1[2] = @elapsed for j = 1:(csc.n)
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
    times1[3] = @elapsed begin
	    colptr[csc.n + 1] = inz
	   	# Julia 1.7 wants this correct
		resize!(rowval, inz - 1)
		resize!(nzval, inz - 1)
		SparseMatrixCSC{Tv, Ti}(csc.m, csc.n, colptr, rowval, nzval)
	end	
    
    return times1
end


"""
`function pe_p_time(lnk::SparseMatrixLNK{Tv, Ti}, csc::SparseMatrixCSC, nt, nt2) where {Tv, Ti <: Integer}`

Timed version of `pe_p`.
"""
function pe_p_time(lnk::SparseMatrixLNK{Tv, Ti},
                 csc::SparseMatrixCSC, nt, nt2) where {Tv, Ti <: Integer}
    @assert(csc.m==lnk.m)
    @assert(csc.n==lnk.n)
    
    times1 = zeros(6+nt)

	times1[1] = @elapsed begin
		nc_nt = Int(csc.n/nt)

		# Detect the maximum column length of lnk
		lnk_maxcols = zeros(Ti, nt)
		xnnzs  = [csc.colptr[tid*nc_nt+1]-csc.colptr[(tid-1)*nc_nt+1] for tid=1:nt]
	end
	
	#nt2 = 4
	fac = Int(nt/nt2)
	times1[2] = @elapsed @threads :static for ltid=1:nt2
		for tid0=1:fac
			tid = (ltid-1)*fac+tid0
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
	end
	
	
	times1[3] =@elapsed	begin
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
	end
    # loop over all columns
    #@threads :static 
    for tid=1:nt
    	
    	inz = 1
    	times1[3+tid] = @elapsed for j=(tid-1)*nc_nt+1:tid*nc_nt
    		
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
    
    times1[3+nt+1] = @elapsed begin
		j0 = 0
		for tid=2:nt
			j0 += colptrs[tid-1][nc_nt+1]-1
			colptrs[1][(tid-1)*nc_nt+1:tid*nc_nt] = colptrs[tid][1:end-1].+j0
		end
		colptrs[1][end] = j0+colptrs[nt][nc_nt+1]
	end
	
	times1[3+nt+2] = @elapsed begin
		#j0 = 0
		@threads :static for tid=2:nt
			j1 = colptrs[1][(tid-1)*nc_nt+1] #colptrs[tid-1][nc_nt+1]-1
			j2 = j1+colptrs[tid][end]-2 # colptrs[1][tid*nc_nt+1]-1
			rowvals[1][j1:j2] = view(rowvals[tid], 1:colptrs[tid][end]-1) #rowvals[tid][1:colptrs[tid][end]-1]
			nzvals[1][j1:j2]  = view(nzvals[tid], 1:colptrs[tid][end]-1) #nzvals[tid][1:colptrs[tid][end]-1]
		end
		
	end
    
    times1[3+nt+3] = @elapsed begin
		resize!(rowvals[1], colptrs[1][end] - 1)
		resize!(nzvals[1], colptrs[1][end] - 1)
		SparseMatrixCSC{Tv, Ti}(csc.m, csc.n, colptrs[1], rowvals[1], nzvals[1])
	end
	return [times1[1], times1[2], times1[3], times1[4:3+nt], times1[end-2], times1[end-1], times1[end]] #times1
end


"""
`function pe_p2_time(lnk::SparseMatrixLNK{Tv, Ti}, csc::SparseMatrixCSC, nt)::SparseMatrixCSC where {Tv, Ti <: Integer}`

Timed version `pe_p2`.
"""
function pe_p2_time(lnk::SparseMatrixLNK{Tv, Ti},
                 csc::SparseMatrixCSC, nt)::SparseMatrixCSC where {Tv, Ti <: Integer}
    @assert(csc.m==lnk.m)
    @assert(csc.n==lnk.n)
    #@info csc.n
    
    m = 3
    times = zeros(m*nt)
    
    stat = 0
    readys = zeros(Ti, nt)
    
    gxnnz = nnz(csc) + nnz(lnk)
    gcolptr = Vector{Ti}(undef, csc.n + 1)
    gcolptr[1] = 1
    growval = Vector{Ti}(undef, gxnnz)
    gnzval = Vector{Tv}(undef, gxnnz)
    nc_nt = Int(csc.n/nt)
    
    in_csc_col(jcsc, j) = (nnz(csc) > zero(Ti)) && (jcsc < csc.colptr[j + 1])
	in_lnk_col(jlnk, l_lnk_col) = (jlnk <= l_lnk_col)
	
	@threads :static for tid=1:nt
		times[tid] = @elapsed begin
			xnnz  = csc.colptr[tid*nc_nt+1]-csc.colptr[(tid-1)*nc_nt+1]
			
			lnk_maxcol = zero(Ti)
			
			inz = 1
			for j=(tid-1)*nc_nt+1:tid*nc_nt
				lcol = zero(Ti)
				k = j
				while k > 0
					lcol += 1
				   	k = lnk.colptr[k]
			   	end
			   	xnnz += lcol
			   	lnk_maxcol = max(lcol, lnk_maxcol)
			end
		
			#colptr = Vector{Ti}(undef, nc_nt + 1) for tid=1:nt]
			rowval = Vector{Ti}(undef, xnnz) #for tid=1:nt]
			nzval  = Vector{Tv}(undef, xnnz) #for tid=1:nt]

			# pre-allocate column  data
			col = [ColEntry{Tv, Ti}(0, zero(Tv)) for i = 1:lnk_maxcol] # for tid=1:nt]

			#inzs = ones(Ti, nt) # counts the nonzero entries in the new matrix
			inz = 1
		end	
		times[tid+nt] = @elapsed for j=(tid-1)*nc_nt+1:tid*nc_nt
    		
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

		    if j != (tid-1)*nc_nt+1
			    gcolptr[j] = inz
			end
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
    	
    	times[tid+2*nt] = @elapsed begin
			gcolptr[tid*nc_nt+1] = inz  # s#[tid][nc_nt + 1] = inz
			readys[tid] = 1
			while (stat != tid-1) || (sum(readys) != nt)
			end
			s
			j1 = gcolptr[(tid-1)*nc_nt + 1]
			j2 = j1+inz-2  #gcolptr[tid*nc_nt+1]-1
			
			
			if tid == 1
			else
				gcolptr[(tid-1)*nc_nt+2:end] .+= gcolptr[(tid-1)*nc_nt+1]-gcolptr[(tid-2)*nc_nt+1]
			end
			
			stat = tid
			
			growval[j1:j2] = view(rowval, 1:j2-j1+1) 
			gnzval[j1:j2]  = view(nzval, 1:j2-j1+1) 
			
		end
    	
    end
    
    @info times[1:nt]
    @info times[nt+1:2*nt]
    @info (times[1:nt]+times[nt+1:2*nt])
    @info times[2*nt+1:3*nt]
    @info (times[1:nt]+times[nt+1:2*nt]+times[2*nt+1:3*nt])
    @info maximum(times[1:nt]), maximum(times[nt+1:2*nt]), maximum(times[2*nt+1:3*nt])
    @info minimum(times[1:nt]), minimum(times[nt+1:2*nt]), minimum(times[2*nt+1:3*nt])
    
    
    #colptr[csc.n + 1] = inz
    # Julia 1.7 wants this correct
    resize!(growval, gcolptr[end] - 1)
    resize!(gnzval, gcolptr[end] - 1)
    SparseMatrixCSC{Tv, Ti}(csc.m, csc.n, gcolptr, growval, gnzval)
end


### compare and benchmark functions



"""
`function compare_seq_to_parallel(num, n, k, nt2)`

1. Use the sequential `pe_s` and the parallel `pe_p` functions for `CSC + LNK` and compare the results.
2. Compare the time and allocations they take / make.
3. `pe_p2` is not included yet!

Random matrices are used. The number of threads used is the number Julia is started with (e.g. `julia -t 4`). `n` has to be an integer multiple of the number of threads. \\
num: Samplesize for the benchmark. \\
n: Both the LNK and the CSC matrix are square n x n matrices. \\
k: Number of changes in the matrix, k <= number of non-zero entries \\
nt2: Number of threads for smaller loops in `pe_p`.
"""
function compare_detail(num, n, k, nt2)
	B, C = create_matrix2(n, k)
	C1 = pe_s(B, C)
	C2 = pe_p(B, C, nt, nt2)
	C3 = pe_p2(B, C, nt)
	
	
	

	println("nt2 = $nt2")
	println("Length of (CSC_sequential - CSC_parallel).nzval = $(length((C1-C2).nzval))")
	println("Length of (CSC_sequential - CSC_parallel2).nzval = $(length((C1-C3).nzval))")
	
	timeso1 = zeros(num, 3)
	timeso2 = zeros(num, 7)
	timeso3 = zeros(num)
	for i=1:num
		@info i
		B, C = create_matrix2(n, k)
		timeso1[i,:] = pe_s_time(B, C)
		tmp = pe_p_time(B, C, nt, nt2)
		timeso2[i,1] = tmp[1] #pe_p_time(B, C, nt)
		timeso2[i,2] = tmp[2] #pe_p_time(B, C, nt)
		timeso2[i,3] = tmp[3] #pe_p_time(B, C, nt)
		timeso2[i,4] = minimum(tmp[4])
		timeso2[i,5] = tmp[5]
		timeso2[i,6] = tmp[6]
		timeso2[i,7] = tmp[7]
		timeso3[i] = maximum(tmp[4])
		
		
		GC.gc()
	end
	
	
	println("pre ", minimum(timeso1[:,1]), ", loop ", minimum(timeso1[:,2]), ", post ", minimum(timeso1[:,3]), " loop/nt ", (minimum(timeso1[:,2])/nt))
	println("pre1 ", minimum(timeso2[:,1]), "pre2 ", minimum(timeso2[:,2]), "pre3 ", minimum(timeso2[:,3]), ", loop ", minimum(timeso2[:,4]), ", colptr ", minimum(timeso2[:,5]), ", vals ", minimum(timeso2[:,6]), ", post ", minimum(timeso2[:,7]), " short ", minimum(timeso3))
	
	
	times3 = zeros(num)
	times4 = zeros(num)
	times5 = zeros(num)
	for i=1:num
		B, C = create_matrix2(n, k)
		times3[i] = @elapsed pe_s(B, C)
		times4[i] = @elapsed pe_p(B, C, nt, nt2)
		times5[i] = @elapsed pe_p2(B, C, nt)
		@warn i, times3[i], times4[i], times5[i]
		
		GC.gc()
	end
	
	
	B, C = create_matrix2(n, k)
	a1 = @allocated pe_s(B, C)
	a2 = @allocated pe_p(B, C, nt, nt2)
	a3 = @allocated pe_p2(B, C, nt)
	#timeso2_short = shortening(timeso3, nt)
	
	
	
	println("Time  for sequential LNK+CSC: $(minimum(times3)) seconds")
	println("Time  for parallel   LNK+CSC: $(minimum(times4)) seconds")
	println("Time  for parallel2  LNK+CSC: $(minimum(times5)) seconds")
	println("Allocs in sequential LNK+CSC: $a1 bytes")
	println("Allocs in parallel   LNK+CSC: $a2 bytes")
	println("Allocs in parallel2  LNK+CSC: $a3 bytes")
	
end


"""
`function compare_quick_all(num, n, k, nt2)`

1. Use the sequential `pe_s` and the parallel `pe_p` & `pe_p2` functions for `CSC + LNK` and compare the results.
2. Compare the time and allocations they take / make.
3. `pe_p2` is not included yet!

Random matrices are used. The number of threads used is the number Julia is started with (e.g. `julia -t 4`). `n` has to be an integer multiple of the number of threads. \\
num: Samplesize for the benchmark. \\
n: Both the LNK and the CSC matrix are square n x n matrices. \\
k: Number of changes in the matrix, k <= number of non-zero entries \\
nt2: Number of threads for smaller loops in `pe_p`.
"""
function compare_quick_all(num, n, k, nt2)
	B, C = create_matrix2(n, k)
	C1 = pe_s(B, C)
	C2 = pe_p(B, C, nt, nt2)
	C3 = pe_p2(B, C, nt)
	
	println("Length of (CSC_sequential - CSC_parallel).nzval  = $(length((C1-C2).nzval))")
	println("Length of (CSC_sequential - CSC_parallel2).nzval = $(length((C1-C3).nzval))")
	times3 = zeros(num)
	times4 = zeros(num)
	times5 = zeros(num)
	for i=1:num
		B, C = create_matrix2(n, k)
		times3[i] = @elapsed pe_s(B, C)
		times4[i] = @elapsed pe_p(B, C, nt, nt2)
		times5[i] = @elapsed pe_p2(B, C, nt)
		@info i, times3[i], times4[i], times5[i]
		
		GC.gc()
	end
	
	B, C = create_matrix2(n, k)
	a1 = @allocated pe_s(B, C)
	a2 = @allocated pe_p(B, C, nt, nt2)
	a3 = @allocated pe_p2(B, C, nt)
	#timeso2_short = shortening(timeso3, nt)
	
	
	
	println("Time for sequential  LNK+CSC: $(minimum(times3)) seconds")
	println("Time  for parallel   LNK+CSC: $(minimum(times4)) seconds")
	println("Time  for parallel2  LNK+CSC: $(minimum(times5)) seconds")
	println("Allocs in sequential LNK+CSC: $a1 bytes")
	println("Allocs in parallel   LNK+CSC: $a2 bytes")
	println("Allocs in parallel2  LNK+CSC: $a3 bytes")
	
end

"""
`function test_pe_s(num, n, k)`

Runs `pe_s` num times with n x n matrices with k nonzero entries.
"""
function test_pe_s(num, n, k)
	B, C = create_matrix2(n, k)
	C1 = pe_s(B, C)
	
	#println("Length of (CSC_sequential - CSC_parallel2).nzval = $(length((C1-C3).nzval))")
	times = zeros(num)
	for i=1:num
		B, C = create_matrix2(n, k)
		@info i
		times[i] = @elapsed pe_s(B, C)
		@info times[i]
		GC.gc()
	end
	
	println(minimum(times), ", ", (sum(times)/num), ", ", maximum(times))	
end

"""
`function test_pe_p2(num, n, k)`

Runs `pe_p2` num times with n x n matrices with k nonzero entries.
"""
function test_pe_p2(num, n, k)
	B, C = create_matrix2(n, k)
	#C1 = pe_s(B, C)
	C3 = pe_p2(B, C, nt)
	
	#println("Length of (CSC_sequential - CSC_parallel2).nzval = $(length((C1-C3).nzval))")
	times = zeros(num)
	for i=1:num
		B, C = create_matrix2(n, k)
		@info i
		times[i] = @elapsed pe_p2(B, C, nt)
		@info times[i]
		GC.gc()
	end
	
	println(minimum(times), ", ", (sum(times)/num), ", ", maximum(times))	
end


"""
`function create_matrix(n, k)`

SLOW!, use `create_matrix2` instead!
Creates two random n x n matrices with k nonzero entries, one matrix is in LNK and one in CSC format.
"""	
function create_matrix(n, k)
	A = SparseMatrixLNK{Float64, Int32}(n, n)
	B = SparseMatrixLNK{Float64, Int32}(n, n)

	for i=1:k
		A[rand(1:n), rand(1:n)] += 1.0
		B[rand(1:n), rand(1:n)] += 1.0
	end

	C = SparseMatrixCSC(A)
	B, C
end	

"""
`function create_matrix2(n, k)`

Creates two random n x n matrices with k nonzero entries, one matrix is in LNK and one in CSC format.
This function uses `sprand_sdd!`.
"""
function create_matrix2(n, k)
	nnzrow = Int(k/n)
	A = SparseMatrixLNK{Float64, Int32}(n, n)
	B = SparseMatrixLNK{Float64, Int32}(n, n)
	sprand_sdd!(A; nnzrow=nnzrow)
	sprand_sdd!(B; nnzrow=nnzrow)
	C = SparseMatrixCSC(A)
	B, C
end	

nt2 = empirical_nt2()


test_pe_s(5, 1_000_000, 10_000_000)
test_pe_p2(5, 1_000_000, 10_000_000)
compare_quick_all(5, 1_000_000, 10_000_000, nt2)
