using SparseArrays
using ExtendableSparse
using ExtendableGrids
using Metis
using Base.Threads

nt = nthreads()

include(path*"preparatory.jl")
include(path*"assembly.jl")
include(path*"csc_assembly.jl")
include(path*"conversion.jl")
include(path*"validation.jl")

validate((10*nt, 10*nt+1))

fct_names = [["CP", "Rm1", "Rm2"], ["CP", "Rm1", "Rm2", "Rm1_dz", "Rm2_dz"], ["CP", "Rm1", "Rm2"]]

endings = ["assembly", "conversion", "cscassembly"]

nms_    = []
grid_   = []
gridm1_ = []
nntsm1_ = []
sm1_    = []
nrm1_   = []
cfpm1_  = []
gridm2_ = []
nntsm2_ = []
sm2_    = []
nrm2_   = []
cfpm2_  = []

As_     = []
Am1_    = []
Am2_    = []
Cs_     = []
Cm1_    = []
Cm2_    = []

"""
`function add_a_grid(nm)`

Add a 2d or 3d grid with (nm)=n,m or (nm)=n,m,l nodes.
Create the preparatory stuff for all assembly functions.
"""
function add_a_grid(nm)
	push!(nms_, nm)
	grid = getgrid(nm)
	push!(grid_, grid)
	push!(As_, da_LNK_cp_sz(grid, nt))
	A1 = copy(As_[end])
	push!(Cs_, CSC_LNKs_s!(A1))

	grid, nnts, s, nr, cfp = preparatory_multi_ps(nm, nt, 1)
	push!(gridm1_, grid)
	push!(nntsm1_, nnts)
	push!(sm1_, s)
	push!(nrm1_, nr)
	push!(cfpm1_, cfp)
	push!(Am1_, da_RLNK_oc_ps_sz(grid, nnts, s, cfp, nt, 1))
	push!(Cm1_, CSC_RLNK_si_oc_ps_dz(Am1_[end], nr, s, nt, 1))
	
	grid, nnts, s, nr, cfp = preparatory_multi_ps(nm, nt, 2)
	push!(gridm2_, grid)
	push!(nntsm2_, nnts)
	push!(sm2_, s)
	push!(nrm2_, nr)
	push!(cfpm2_, cfp)
	push!(Am2_, da_RLNK_oc_ps_sz(grid, nnts, s, cfp, nt, 2))
	push!(Cm2_, CSC_RLNK_si_oc_ps_dz(Am2_[end], nr, s, nt, 2))
	
	#check if file exists already, if not, create it
	for (i,ending) in enumerate(endings)
		fn = pre_name*"_"*tts(nm, nt)*"_"*ending*".txt"
		
		io = open(fn, "a")
		close(io)
		
		io = open(fn, "r")
		if length(readlines(io)) == 0
			#@warn "file does not exist yet"
			close(io)
			io = open(fn, "w")
			write(io, "ss 0\n")
			#write(io, array_to_string(fct_names))
			for fct_name in fct_names[i]
				write(io, fct_name*"\n")
			end		
		end
		
		close(io)
	end
		
end


"""
`function strarray_to_string(strings)`

`["a", "b", "c"] -> "a b c"`
"""
function strarray_to_string(strings) 
	s = strings[1]
	for string in strings[2:end]
		s = s*" "*string
	end
	s
end

"""
`function numarray_to_strarray(nums)`

`[1, 3, 4] -> ["1", "3", "4"]`
"""
function numarray_to_strarray(nums)
	strings = v = Vector{String}(undef, length(nums))
	for (i,x) in enumerate(nums)
		strings[i] = "$x"
	end
	strings
end


"""
`function tts(xx, y)`

Turns `(a,b),y -> "2_str(y)_str(a)_str(b)"` or `(a,b,c),y -> "3_str(y)_str(a)_str(b)_str(c)"`.
"""
function tts(xx, y) #tuple to string
	s = "$(length(xx))_$(y)"
	for x in xx
		s = s*"_$x"
	end
	s
end


### Benchmark assembly functions:
###----------------------------------------------------------------------------------------------

"""
`function bm_da_CP(flag, ss; doallocs=false)`

Assembles the flag'th grid using `cheap parallelization`.
Does this `ss` times and returns the recorded times.
The result can be used for comparison.
"""
function bm_da_CP(flag, ss; doallocs=false)
	grid = grid_[flag]
	t  = zeros(ss)
	for j=1:ss
		t[j] = @elapsed da_LNK_cp_sz(grid, nt)
		GC.gc()
	end
	t
end

function al_da_CP(flag)
	grid = grid_[flag]
	@allocated da_LNK_cp_sz(grid, nt)
end



"""
`function bm_da_Rocm1(flag, ss; doallocs=false)`

Assembles the flag'th grid using `RLNK_ocm1` (overlapping columns using the function for arbitrary levels of separator partition, here only one level: the separator is not partitioned further: nt+1 matrices).
Does this `ss` times and returns the recorded times.
"""
function bm_da_Rm1(flag, ss; doallocs=false)
	gridm1, nntsm1, sm1, nrm1, cfpm1 = gridm1_[flag], nntsm1_[flag], sm1_[flag], nrm1_[flag], cfpm1_[flag]
	depth = 1
	t = zeros(ss)
	for j=1:ss
		t[j] = @elapsed da_RLNK_oc_ps_sz(gridm1, nntsm1, sm1, cfpm1, nt, depth)
		GC.gc()
	end
	t
end

function al_da_Rm1(flag)
	depth = 1
	gridm1, nntsm1, sm1, nrm1, cfpm1 = gridm1_[flag], nntsm1_[flag], sm1_[flag], nrm1_[flag], cfpm1_[flag]
	@allocated da_RLNK_oc_ps_sz(gridm1, nntsm1, sm1, cfpm1, nt, depth)
end


"""
`function bm_da_Rocm2(flag, ss; doallocs=false)`

Assembles the flag'th grid using `RLNK_ocm2` (overlapping columns using the function for arbitrary levels of separator partition, here two levels: 2*nt+1 matrices).
Does this `ss` times and returns the recorded times.
"""
function bm_da_Rm2(flag, ss; doallocs=false)
	gridm2, nntsm2, sm2, nrm2, cfpm2 = gridm2_[flag], nntsm2_[flag], sm2_[flag], nrm2_[flag], cfpm2_[flag]
	depth = 2
	t = zeros(ss)
	for j=1:ss
		t[j] = @elapsed da_RLNK_oc_ps_sz(gridm2, nntsm2, sm2, cfpm2, nt, depth)
		GC.gc()
	end
	t
end	

function al_da_Rm2(flag)
	gridm2, nntsm2, sm2, nrm2, cfpm2 = gridm2_[flag], nntsm2_[flag], sm2_[flag], nrm2_[flag], cfpm2_[flag]
	depth = 2
	@allocated da_RLNK_oc_ps_sz(gridm2, nntsm2, sm2, cfpm2, nt, depth)
end	


### Benchmark conversion functions:
###----------------------------------------------------------------------------------------------

"""
`function bm_c_LNKs_s(ss, grid)`

Benchmarks the conversion of multiple LNK matrices to one CSC matrix
"""
function bm_c_CP(flag, ss)
	grid = grid_[flag]
	
	t = zeros(ss)
	for j=1:ss
		CP = copy(As_[flag])
		t[j] = @elapsed CSC_LNKs_s!(CP)
		GC.gc()
	end
	t
end

function al_c_CP(flag)
	grid = grid_[flag]
	CP = copy(As_[flag])
	@allocated CSC_LNKs_s!(CP)
end

"""
`function bm_c_RLNK_ps1(flag, ss)`

Benchmarks the conversion of an RLNK matrix with depth 1 (i.e. multiple LNK sub-matrices) to a CSC matrix.
"""
function bm_c_Rm1(flag, ss)
	depth = 1
	gridm1, nntsm1, sm1, nrm1, cfpm1 = gridm1_[flag], nntsm1_[flag], sm1_[flag], nrm1_[flag], cfpm1_[flag]
	Roc = Am1_[flag] #da_RLNK_oc_ps(gridm1, nntsm1, sm1, cfpm1, nt, depth)

	t = zeros(ss)
	for j=1:ss
		t[j] = @elapsed CSC_RLNK_si_oc_ps(Roc, nrm1, sm1, nt, depth)
		GC.gc()
	end
	t
end

function al_c_Rm1(flag)
	depth = 1
	gridm1, nntsm1, sm1, nrm1, cfpm1 = gridm1_[flag], nntsm1_[flag], sm1_[flag], nrm1_[flag], cfpm1_[flag]
	Roc = Am1_[flag] #da_RLNK_oc_ps(gridm1, nntsm1, sm1, cfpm1, nt, depth)
	@allocated CSC_RLNK_si_oc_ps(Roc, nrm1, sm1, nt, depth)
end


"""
`function bm_c_RLNK_ps1(flag, ss)`

Benchmarks the conversion of an RLNK matrix with depth 1 (i.e. multiple LNK sub-matrices) to a CSC matrix.
"""
function bm_c_Rm2(flag, ss)
	depth = 2
	gridm2, nntsm2, sm2, nrm2, cfpm2 = gridm2_[flag], nntsm2_[flag], sm2_[flag], nrm2_[flag], cfpm2_[flag]
	Roc = Am2_[flag] #da_RLNK_oc_ps(gridm2, nntsm2, sm2, cfpm2, nt, depth)

	t = zeros(ss)
	for j=1:ss
		t[j] = @elapsed CSC_RLNK_si_oc_ps(Roc, nrm2, sm2, nt, depth)
		GC.gc()
	end
	t
end

function al_c_Rm2(flag)
	depth = 2
	gridm2, nntsm2, sm2, nrm2, cfpm2 = gridm2_[flag], nntsm2_[flag], sm2_[flag], nrm2_[flag], cfpm2_[flag]
	Roc = Am2_[flag] #da_RLNK_oc_ps(gridm2, nntsm2, sm2, cfpm2, nt, depth)
	@allocated CSC_RLNK_si_oc_ps(Roc, nrm2, sm2, nt, depth)
end


"""
`function bm_c_RLNK_ps1_dz(flag, ss)`

Benchmarks the conversion of an RLNK matrix with depth 1 (i.e. multiple LNK sub-matrices) to a CSC matrix.
"""
function bm_c_Rm1_dz(flag, ss)
	depth = 1
	gridm1, nntsm1, sm1, nrm1, cfpm1 = gridm1_[flag], nntsm1_[flag], sm1_[flag], nrm1_[flag], cfpm1_[flag]
	Roc = Am1_[flag] #da_RLNK_oc_ps(gridm1, nntsm1, sm1, cfpm1, nt, depth)

	t = zeros(ss)
	for j=1:ss
		t[j] = @elapsed CSC_RLNK_si_oc_ps_dz(Roc, nrm1, sm1, nt, depth)
		GC.gc()
	end
	t
end

function al_c_Rm1_dz(flag)
	depth = 1
	gridm1, nntsm1, sm1, nrm1, cfpm1 = gridm1_[flag], nntsm1_[flag], sm1_[flag], nrm1_[flag], cfpm1_[flag]
	Roc = Am1_[flag] #da_RLNK_oc_ps(gridm1, nntsm1, sm1, cfpm1, nt, depth)
	@allocated CSC_RLNK_si_oc_ps_dz(Roc, nrm1, sm1, nt, depth)
end


"""
`function bm_c_RLNK_ps1_dz(flag, ss)`

Benchmarks the conversion of an RLNK matrix with depth 1 (i.e. multiple LNK sub-matrices) to a CSC matrix.
"""
function bm_c_Rm2_dz(flag, ss)
	depth = 2
	gridm2, nntsm2, sm2, nrm2, cfpm2 = gridm2_[flag], nntsm2_[flag], sm2_[flag], nrm2_[flag], cfpm2_[flag]
	Roc = Am2_[flag] #da_RLNK_oc_ps(gridm2, nntsm2, sm2, cfpm2, nt, depth)

	t = zeros(ss)
	for j=1:ss
		t[j] = @elapsed CSC_RLNK_si_oc_ps_dz(Roc, nrm2, sm2, nt, depth)
		GC.gc()
	end
	t
end

function al_c_Rm2_dz(flag)
	depth = 2
	gridm2, nntsm2, sm2, nrm2, cfpm2 = gridm2_[flag], nntsm2_[flag], sm2_[flag], nrm2_[flag], cfpm2_[flag]
	Roc = Am2_[flag] #da_RLNK_oc_ps(gridm2, nntsm2, sm2, cfpm2, nt, depth)
	@allocated CSC_RLNK_si_oc_ps_dz(Roc, nrm2, sm2, nt, depth)
end


###	Benchmark CSC assembly functions
### --------------------------------------------------------------------


### Benchmark CSC assembly

"""
`function bm_ca_CP(flag, ss; offset=1)`

"""
function bm_ca_CP(flag, ss; offset=1)
	grid = grid_[flag]
	A = da_LNK_cp_sz(grid, nt)
	C0 = dropzeros(CSC_LNKs_s!(A))
	
	t = zeros(ss)
	for j=1:ss
		C = copy(C0)
		t[j] = @elapsed da_csc_LNK_cp!(C, grid; offset=offset)
		GC.gc()
	end
	t
end

function al_ca_CP(flag; offset=1)
	grid = grid_[flag]
	A = da_LNK_cp_sz(grid, nt)
	C = dropzeros(CSC_LNKs_s!(A))
	@allocated da_csc_LNK_cp!(C, grid; offset=offset)
end

"""
`function bm_ca_Rm1(flag, ss; offset=1)`

"""
function bm_ca_Rm1(flag, ss; offset=1)
	grid, nnts, s, nr, cfp = gridm1_[flag], nntsm1_[flag], sm1_[flag], nrm1_[flag], cfpm1_[flag]
	
	A = Am1_[flag]#da_RLNK_oc_ps(grid, nnts, s, cfp, nt, 1)
	C0 = Cm1_[flag]#dropzeros(CSC_RLNK_si_oc_ps(A, nr, s, nt, 1))
	
	t = zeros(ss)
	for j=1:ss
		C = copy(C0)
		t[j] = @elapsed da_csc_RLNK_oc_ps_sz!(C, grid, nnts, s, cfp, nr, nt, 1; offset=offset)
		GC.gc()
	end
	t
end

function al_ca_Rm1(flag; offset=1)
	grid, nnts, s, nr, cfp = gridm1_[flag], nntsm1_[flag], sm1_[flag], nrm1_[flag], cfpm1_[flag]
	
	A = Am1_[flag]#da_RLNK_oc_ps(grid, nnts, s, cfp, nt, 1)
	C0 = copy(Cm1_[flag]) #dropzeros(CSC_RLNK_si_oc_ps(A, nr, s, nt, 1))
	
	@allocated da_csc_RLNK_oc_ps_sz!(C0, grid, nnts, s, cfp, nr, nt, 1; offset=offset)
	
end

"""
`function bm_ca_Rm2(flag, ss; offset=1)`

"""
function bm_ca_Rm2(flag, ss; offset=1)
	grid, nnts, s, nr, cfp = gridm2_[flag], nntsm2_[flag], sm2_[flag], nrm2_[flag], cfpm2_[flag]
	
	A = Am2_[flag]
	C0 = Cm2_[flag]
	
	t = zeros(ss)
	for j=1:ss
		C = copy(C0)
		t[j] = @elapsed da_csc_RLNK_oc_ps_sz!(C, grid, nnts, s, cfp, nr, nt, 2; offset=offset)
		GC.gc()
	end
	t
end

function al_ca_Rm2(flag; offset=1)
	grid, nnts, s, nr, cfp = gridm2_[flag], nntsm2_[flag], sm2_[flag], nrm2_[flag], cfpm2_[flag]
	
	A = Am2_[flag]
	C0 = copy(Cm2_[flag])
	
	@allocated da_csc_RLNK_oc_ps_sz!(C0, grid, nnts, s, cfp, nr, nt, 2; offset=offset)
	
end



### Run things:
### ------------------------------------------------------------------------

bm_fcts = [[bm_da_CP, bm_da_Rm1, bm_da_Rm2], [bm_c_CP, bm_c_Rm1, bm_c_Rm2, bm_c_Rm1_dz, bm_c_Rm2_dz], [bm_ca_CP, bm_ca_Rm1, bm_ca_Rm2]]
al_fcts = [[al_da_CP, al_da_Rm1, al_da_Rm2], [al_c_CP, al_c_Rm1, al_c_Rm2, al_c_Rm1_dz, al_c_Rm2_dz], [al_ca_CP, al_ca_Rm1, al_ca_Rm2]]

"""
`function run_auto(wo, ss)`

wo: which one should be benchmarked?
ss: sample size
"""
function run_auto(wo, ss)
	for flag=1:length(nms_)
		str_list = []
		fn = pre_name*"_"*tts(nms_[flag], nt)*"_"*endings[wo]*".txt"
		io = open(fn, "r")
	
		for line in readlines(io)
			push!(str_list, line)
			
		end
		
		fcts = bm_fcts[wo] #[bm_ca_CP, bm_ca_Rm1, bm_ca_Rm2]
		
		close(io)
		
		#fn = pre_name*"_"*tts(nms_[flag])*"_assembly.txt"
		io = open(fn, "w")
	
		for (ctr,str) in enumerate(str_list)
			array = split(str, " ")
			if array[1] == "ss"
				total_ss = parse(Int64, array[2])+ss
				write(io, "ss "*"$(total_ss)\n")
			else
				times = fcts[ctr-1](flag, ss)
				array = vcat(array, numarray_to_strarray(times))
				write(io, strarray_to_string(array)*"\n")
				
			end
		end
		
		
		close(io)
	end
	
end


function run_auto_nowrite(;ss=2)
	for wo=1:3
		for flag=1:length(nms_)
			fcts = bm_fcts[wo]
			for f in fcts
				f(flag, ss)
			end
		end
	end	
end



function run_auto_allocs()
	for flag=1:length(nms_)
		str_list = []
		fn = pre_name*"_"*tts(nms_[flag], nt)*"_allocs.txt"
		io = open(fn, "w")
	
		for (i,ending) in enumerate(endings)
			fcts = al_fcts[i]
			write(io, ending*"\n")
			for (j,f) in enumerate(fcts)
				f(flag)
				write(io, fct_names[i][j]*" $(f(flag))\n")
			end
		end
		
		
		close(io)
	end
	
end

