function get_nnnts_and_sortednodesperthread_and_noderegs_from_cellregs_max(cellregs, allcells, start, nn, Ti)
	nt = maximum(cellregs)-1
	#@warn "nt = $nt"

	#loop over each node, get the cellregion of the cell (the one not in the separator) write the position of that node inside the cellregions sorted ranking into a long vector
	nnts = zeros(Ti, nt+1)
	nnts_max = zeros(Ti, nt+1)
	noderegs = zeros(Ti, nn)
	sortednodesperthread = zeros(Ti, nn)
	sortednodesperthread_max = zeros(Ti, nn)
	noderegs_max_tmp = 0
	for j=1:nn
		cells = @view allcells[start[j]:start[j+1]-1]
		noderegs[j] = minimum(cellregs[cells]) #minimum, s.t the separator is only choosen if the node is just in the separator
		noderegs_max_tmp = noderegs[j]
		nnts[noderegs[j]] += 1
		sortednodesperthread[j] = nnts[noderegs[j]]
		if noderegs[j] != maximum(cellregs[cells])
			noderegs[j] *= -1
			noderegs_max_tmp = nt+1
		end
		nnts_max[noderegs_max_tmp]  += 1
		sortednodesperthread_max[j] = nnts_max[noderegs_max_tmp]
		#end
	end
	
	nnts, sortednodesperthread, noderegs, nnts_max, sortednodesperthread_max
end

function get_nnnts_and_sortednodesperthread_and_noderegs_from_cellregs_ps(cellregs, allcells, start, nn, Ti, nt)
		
	num_matrices = maximum(cellregs)
	depth = Int(floor((num_matrices-1)/nt))

	#loop over each node, get the cellregion of the cell (the one not in the separator) write the position of that node inside the cellregions sorted ranking into a long vector
	#nnts = [zeros(Ti, nt+1) for i=1:depth+1]
	nnts = zeros(Ti, depth*nt+1)
	sortednodesperthread = zeros(Ti, (depth+1, nn)) #[zeros(Ti, nn) for i=1:depth]
	#noderegs_max_tmp = 0
	noderegions = zeros(Ti, (depth+1, nn))
	

	for j=1:nn
		cells = @view allcells[start[j]:start[j+1]-1]
		sortedcellregs = unique(sort(cellregs[cells])) #.unique
		#nr1 = minimum(cellregs[cells])
		#nr2 = maximum(cellregs[cells])

		#i1 = Int(floor(nr1/nt))
		#i2 = Int(floor(nr2/nt))

		for cr in sortedcellregs #i=i1:i2
			i = Int(ceil(cr/nt))
			#nnts[i][cr] += 1
			nnts[cr] += 1
			sortednodesperthread[i,j] = nnts[cr] #nnts[i][cr]
			noderegions[i,j] = cr
		end

	end
	nnts, sortednodesperthread, noderegions
end

function separate!(cellregs, nc, ACSC, nt, level0, ctr_sepanodes)
	sepanodes = findall(x->x==nt+1, cellregs)

	indptr = collect(1:nc+1)
	indices = zeros(Int64, nc)
	rowval = zeros(Int64, nc)

	indptrT = collect(1:ctr_sepanodes+1)
	indicesT = zeros(Int64, ctr_sepanodes)
	rowvalT = zeros(Int64, ctr_sepanodes)

	for (i,j) in enumerate(sepanodes)
		indices[j] = i
		indicesT[i] = j
		rowval[j]  = 1
		rowvalT[i] = 1
	end

	R = SparseMatrixCSC(ctr_sepanodes, nc, indptr, indices, rowval)
	RT = SparseMatrixCSC(nc, ctr_sepanodes, indptrT, indicesT, rowvalT)
	prod = ACSC*dropzeros(RT)
	RART = dropzeros(R)*ACSC*dropzeros(RT)
	
	partition2 = Metis.partition(RART, nt)
	cellregs2 = copy(partition2)

	ctr_sepanodes = 0
	for (i,j) in enumerate(sepanodes)
		rows = RART.rowval[RART.colptr[i]:(RART.colptr[i+1]-1)]
		cellregs[j] = level0*nt + cellregs2[i]
		if minimum(partition2[rows]) != maximum(partition2[rows])
			cellregs[j] = (level0+1)*nt+1
			ctr_sepanodes += 1
		end
	end

	RART, ctr_sepanodes
end

function grid_to_graph!(grid::ExtendableGrid, nt::Int)
	#For each node I want to know in which cells it is contained
	#Such that I can connect each pair of cells from one node with an edge (Remark: in the graph we are building, the grid-cells are the vertices)

	#The matrix assembly and the computation of the separator can easily be parallelized. The sppedup is questionable though.
	
	A = SparseMatrixLNK{Int64, Int64}(num_cells(grid), num_cells(grid))
	#Count how many cells contain the `node_id`-th node 
	number_cells_per_node = zeros(Int64, num_nodes(grid))
	for j=1:num_cells(grid)
		for node_id in grid[CellNodes][:,j]
			number_cells_per_node[node_id] += 1
		end
	end

	#allcells = [cell1_from_node1, cell2_from_node1, ..., cell6_from_node1, cell1_from_node2, ..., cell6_from_node1000]
	allcells = zeros(Int64, sum(number_cells_per_node))
	#The indices at which the transition from node1 to node2 etc in `allcells` occur are stored in `start` (basically cumsum(number_cells_per_node))
	start = ones(Int64, num_nodes(grid)+1)
	start[2:end] += cumsum(number_cells_per_node)
	#counter of cells per node
	number_cells_per_node .= 0
	for j=1:num_cells(grid)
		for node_id in grid[CellNodes][:,j]
			allcells[start[node_id] + number_cells_per_node[node_id]] = j
			number_cells_per_node[node_id] += 1
		end
	end

#Iterate over all nodes, store the cells that contain this node in `cells` and then iterate over all cell pairs to connect them in the adjacency matrix
	for j=1:num_nodes(grid)
		cells = @view allcells[start[j]:start[j+1]-1]
		for (i,id1) in enumerate(cells)
			for id2 in cells[i+1:end]
				A[id1,id2] = 1
				A[id2,id1] = 1
			end
		end
		
	end

	#The LNK adjacency matrix is converted to CSC
	ACSC = SparseArrays.SparseMatrixCSC(A)
	#The graph is partitioned, such that the number of partitions is equal to the number of threads `nt`
	
	partition = Metis.partition(ACSC, nt)
	cellregs  = copy(partition)

	

	#The goal is to find all cells (i.e. vertices in the graph) that are connected to cells of different partition-number i.e. cellregion. For each of these cells, we change the partition number to `nt+1`, i.e. they are in the separator
	#We look at all cells the `j`-th cell is connected to, by looking at all non-zero entries in the `j`-th column.
	for j=1:num_cells(grid)
		rows = ACSC.rowval[ACSC.colptr[j]:(ACSC.colptr[j+1]-1)]
		if minimum(partition[rows]) != maximum(partition[rows])
			cellregs[j] = nt+1
		end
	end

	grid[CellRegions] = cellregs
	return allcells, start
end


function grid_to_graph_ps_multi!(grid, nt, depth)
	A = SparseMatrixLNK{Int64, Int64}(num_cells(grid), num_cells(grid))
	number_cells_per_node = zeros(Int64, num_nodes(grid))
	for j=1:num_cells(grid)
		for node_id in grid[CellNodes][:,j]
			number_cells_per_node[node_id] += 1
		end
	end
	allcells = zeros(Int64, sum(number_cells_per_node))
	start = ones(Int64, num_nodes(grid)+1)
	start[2:end] += cumsum(number_cells_per_node)
	number_cells_per_node .= 0
	for j=1:num_cells(grid)
		for node_id in grid[CellNodes][:,j]
			allcells[start[node_id] + number_cells_per_node[node_id]] = j
			number_cells_per_node[node_id] += 1
		end
	end

	for j=1:num_nodes(grid)
		cells = @view allcells[start[j]:start[j+1]-1]
		for (i,id1) in enumerate(cells)
			for id2 in cells[i+1:end]
				A[id1,id2] = 1
				A[id2,id1] = 1
			end
		end	
	end

	ACSC = SparseArrays.SparseMatrixCSC(A)
	
	partition = Metis.partition(ACSC, nt)
	cellregs  = copy(partition)
	
	ctr_sepanodes = 0
	for j=1:num_cells(grid)
		rows = ACSC.rowval[ACSC.colptr[j]:(ACSC.colptr[j+1]-1)]
		if minimum(partition[rows]) != maximum(partition[rows])
			cellregs[j] = nt+1
			ctr_sepanodes += 1
		end
	end
	RART = ACSC
	for level=1:depth-1
		RART, ctr_sepanodes = separate!(cellregs, num_cells(grid), RART, nt, level, ctr_sepanodes)
	end

			
	grid[CellRegions] = cellregs
	
	return allcells, start
end





function vvcons(Ti, lengths)
	x::Vector{Vector{Ti}} = [zeros(Ti, i) for i in lengths]
	return x
end

function bettercellsforpart(xx, upper)
	ctr = zeros(Int64, upper)
	for x in xx
		ctr[x] += 1
	end
	cfp = vvcons(Int64, ctr)
	ctr .= 1
	for (i,x) in enumerate(xx)
		cfp[x][ctr[x]] = i
		ctr[x] += 1
	end
	cfp
end

function getgrid(nm)
	if length(nm) == 2
		n,m = nm
		xx = collect(LinRange(0.0, 1.0, n))
		yy = collect(LinRange(0.0, 1.0, m))
		grid = simplexgrid(xx, yy)
	else 
		n,m,l = nm
		xx = collect(LinRange(0.0, 1.0, n))
		yy = collect(LinRange(0.0, 1.0, m))
		zz = collect(LinRange(0.0, 1.0, l))
		grid = simplexgrid(xx, yy, zz)
	end
	grid
end

function preparatory_max(nm, nt)
	grid = getgrid(nm)
	
	allcells, start = grid_to_graph!(grid, nt)

	nnts, s, nr, nnts2, s2 = get_nnnts_and_sortednodesperthread_and_noderegs_from_cellregs_max(
		grid[CellRegions], allcells, start, num_nodes(grid), Int32
		)
	return grid, nnts, s, nr, nnts2, s2, bettercellsforpart(grid[CellRegions], nt+1)
end

function preparatory_multi_ps(nm, nt, depth)
	grid = getgrid(nm)
	
	t = @elapsed ((allcells, start) = grid_to_graph_ps_multi!(grid, nt, depth))
	#@warn "grid to graph $t"
	t = @elapsed ((nnts, s, nr) = get_nnnts_and_sortednodesperthread_and_noderegs_from_cellregs_ps(
		grid[CellRegions], allcells, start, num_nodes(grid), Int32, nt
	))
	#@warn "get... $t"
	return grid, nnts, s, nr, bettercellsforpart(grid[CellRegions], depth*nt+1)
end

function preparatory_CP(nm, nt)
	grid = getgrid(nm)
	
	return grid, nt
end
