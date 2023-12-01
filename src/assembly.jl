"""
`function fr(x)`

fr stands for fake random and it is supposed to imitate a random number generator, but reproducible.
It just returns `(x%3)-1`. Used in matrix assembly to have different nonzero entries in the original and the 'new' CSC matrix.
"""
function fr(x)
	(x%3)-1
end


"""
`function da_RLNK_oc_ps_sz(grid, nnts, s, cellsforpart, nt, depth)`

Dummy assembly (da) of Reduced LNK matrices with overlapping columns (at thi spoint, this is the standard) with a partitioned separator (if `depth`>1) and some zeros (since the off-diagonal entries have the sign -1, 0 or 1).
"""
function da_RLNK_oc_ps_sz(grid, nnts, s, cellsforpart, nt, depth)
	K = size(grid[CellNodes])[1]
	As = [SparseMatrixLNK{Float64, Int32}(num_nodes(grid), nnts[tid]) for tid=1:depth*nt+1]

	for level=1:depth
		@threads :static for tid=(level-1)*nt+1:level*nt #smth=1:nt
			for icell in cellsforpart[tid] 
				for i=1:K
					inode = grid[CellNodes][i,icell]
					As[tid][inode, s[level,inode]] += 3.0
					for j=i+1:K
						jnode = grid[CellNodes][j,icell]
						v = fr(inode+jnode)*0.5
						As[tid][inode, s[level,jnode]] += v
						As[tid][jnode, s[level,inode]] += v
					end
				end
			end
		end
	end

	for icell in cellsforpart[depth*nt+1] 
		for i=1:K
			inode = grid[CellNodes][i,icell]
			As[depth*nt+1][inode, s[depth+1,inode]] += 3.0
			for j=i+1:K
				jnode = grid[CellNodes][j,icell]
				v = fr(inode+jnode)*0.5
				As[depth*nt+1][inode, s[depth+1,jnode]] += v
				As[depth*nt+1][jnode, s[depth+1,inode]] += v
			end
		end
	end
	
 	As

end

"""
`function da_LNK_cp_sz(grid::ExtendableGrid, nt::Integer)`

Dummy assembly (da) of an LNK using 'using cheap parallelization' with some zeros (since the off-diagonal entries have the sign -1, 0 or 1).
"""
function da_LNK_cp_sz(grid::ExtendableGrid, nt::Integer)
	K = size(grid[CellNodes])[1]
	nc_nt = Int(num_cells(grid)/nt)
	cfp = [collect((tid-1)*nc_nt+1:tid*nc_nt) for tid=1:nt]
	As = [SparseMatrixLNK{Float64, Int32}(num_nodes(grid), num_nodes(grid)) for tid=1:nt]

	@threads :static for tid=1:nt
		for icell in cfp[tid] 
			for i=1:K
				inode = grid[CellNodes][i,icell]
				As[tid][inode, inode] += 3.0
				for j=i+1:K
					jnode = grid[CellNodes][j,icell]
					v = fr(inode+jnode)*0.5	
					As[tid][inode,jnode] += v #(-0.5)
					As[tid][jnode,inode] += v #(-0.5)
				end
			end
		end
	end

	As
end



