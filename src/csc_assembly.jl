include(path*"assembly.jl")

"""
`function entryexists(CSC, i, j)`

OLD: find out if CSC already has an nonzero entry at i,j.
`entryexists2` should be used instead.
"""

function entryexists(CSC, i, j)
	ids = CSC.colptr[j]:(CSC.colptr[j+1]-1)
	i in CSC.rowval[ids]
end

"""
`function entryexists2(CSC, i, j)`

Find out if CSC already has an nonzero entry at i,j without any allocations
"""
function entryexists2(CSC, i, j) #find out if CSC already has an nonzero entry at i,j
	#vals = 
	#ids = CSC.colptr[j]:(CSC.colptr[j+1]-1)
	i in view(CSC.rowval, CSC.colptr[j]:(CSC.colptr[j+1]-1))
end

"""
`function updatentryCSC!(CSC, i, j, v)`

CSC[i,j] += v, without any allocations
"""
function updatentryCSC!(CSC, i, j, v)
	p = CSC.colptr[j]
	while p<CSC.colptr[j+1]
		r = CSC.rowval[p]
		if i==r
			CSC.nzval[p] += v
			return
		end
		p += 1
	end
end

"""
`function da_csc_LNK_s!(C, grid; offset=0)`

Dummy assembly in the CSC matrix sequentially (this is the continuation of the basic LNK assembly)
"""
function da_csc_LNK_s!(C, grid; offset=0)
	C.nzval .= 0
	for icell=1:num_cells(grid)
		for (i,inode) in enumerate(grid[CellNodes][:,icell])
			C[inode,inode] += 3.0
			for jnode in grid[CellNodes][i+1:end,icell]
				v = fr(inode+jnode+offset)*0.5	
				C[inode,jnode] += v #(-0.5)
				C[jnode,inode] += v #(-0.5)
			end
		end
	end
end


"""
`function da_csc_LNK_cp!(C, grid; offset=0)`

Dummy assembly in the CSC matrix cheaply parallelized
"""
function da_csc_LNK_cp!(C, grid; offset=0)
	A_backup = [SparseMatrixLNK{Float64, Int32}(num_nodes(grid), num_nodes(grid)) for i=1:nt]
	nc_nt = Int(num_cells(grid)/nt)
	cfp = [collect((tid-1)*nc_nt+1:tid*nc_nt) for tid=1:nt]
	K = size(grid[CellNodes])[1]
	
	C.nzval .= 0
	
	@threads :static for tid=1:nt
		for icell in cfp[tid] 
			#for (i,inode) in enumerate(grid[CellNodes][:,icell])
			for i=1:K
				inode = grid[CellNodes][i,icell]
				if entryexists2(C, inode, inode)
					updatentryCSC!(C, inode, inode, 3.0)
				else
					A_backup[tid][inode,inode] += 3.0
				end
				
				#for jnode in grid[CellNodes][i+1:end,icell]
				for j=i+1:K
					jnode = grid[CellNodes][j,icell]
					v = fr(inode+jnode+offset)*0.5	
					if entryexists2(C, inode, jnode)
						updatentryCSC!(C, inode, jnode, v)
						updatentryCSC!(C, jnode, inode, v)
					else
						A_backup[tid][inode,jnode] += v
						A_backup[tid][jnode,inode] += v
					end
					
				end
			end
		end
	end

	
	
	for tid=2:nt
		add_LNKs!(A_backup[1], A_backup[tid])
	end
	A_backup[1] + C
	
end

"""
`function fill_dummy_zeros!(RLNKs)`

Fills the first row of the matrices with zeros, if this is not done there might be problems with CSC conversion.
`RLNKs` is a vector of LNK matrices of reduced size.
"""
function fill_dummy_zeros!(RLNKs)
	@threads :static for tid=1:length(RLNKs) #depth*nt+1
		for j=1:RLNKs[tid].n
			RLNKs[tid][1,j] = 1.0
			RLNKs[tid][1,j] -= 1.0
		end
	end
end


"""
`function da_csc_RLNK_oc_ps_sz!(C, grid, nnts, s, cellsforpart, nr, nt, depth; offset=0)`

Dummy assemble in the CSC matrix based on a grid partition with some depth
"""
function da_csc_RLNK_oc_ps_sz!(C, grid, nnts, s, cellsforpart, nr, nt, depth; offset=0)
	K = size(grid[CellNodes])[1] #number of nodes in a cell, obviously only forks for simplexgrids
	C.nzval .= 0
	A_backup = [SparseMatrixLNK{Float64, Int32}(num_nodes(grid), nnts[tid]) for tid=1:depth*nt+1]

	fill_dummy_zeros!(A_backup)

	for level=1:depth
		@threads :static for tid=(level-1)*nt+1:level*nt #smth=1:nt
			for icell in cellsforpart[tid]
				#for (i,inode) in enumerate(grid[CellNodes][:,icell])
				for i=1:K
					inode = grid[CellNodes][i,icell]
					
					if entryexists2(C, inode, inode)
						updatentryCSC!(C, inode, inode, 3.0)
					else
						A_backup[tid][inode, s[level,inode]] += 3.0
					end
					#for jnode in grid[CellNodes][i+1:end,icell]
					for j=i+1:K
						jnode = grid[CellNodes][j,icell]
						v = fr(inode+jnode+offset)*0.5	
						if entryexists2(C, inode, jnode)
							updatentryCSC!(C, inode, jnode, v)
							updatentryCSC!(C, jnode, inode, v)
						else
							A_backup[tid][inode, s[level,jnode]] += v
							A_backup[tid][jnode, s[level,inode]] += v
						end
					end
				end
			end
		end
	end

	for icell in cellsforpart[depth*nt+1] 
		#for (i,inode) in enumerate(grid[CellNodes][:,icell])
		for i=1:K
			inode = grid[CellNodes][i,icell]
			if entryexists2(C, inode, inode)
				updatentryCSC!(C, inode, inode, 3.0)
			else
				A_backup[depth*nt+1][inode, s[depth+1,inode]] += 3.0
			end
			#for jnode in grid[CellNodes][i+1:end,icell]
			for j=i+1:K
				jnode = grid[CellNodes][j,icell]
				v = fr(inode+jnode+offset)*0.5	
				if entryexists2(C, inode, jnode)
					updatentryCSC!(C, inode, jnode, v)
					updatentryCSC!(C, jnode, inode, v)
				else
					A_backup[depth*nt+1][inode, s[depth+1,jnode]] += v
					A_backup[depth*nt+1][jnode, s[depth+1,inode]] += v
				end
			end
		end
	end
	
	C+CSC_RLNK_si_oc_ps_dz(A_backup, nr, s, nt, depth)
	
end
