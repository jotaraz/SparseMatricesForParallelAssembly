"""
`function validate(nm::Tuple{Integer,Integer})`

This used three differnet ways to compute the same matrices of an n x m (or n x m x l) grid.
If it outputs 4x`Float64[]`, everything is successful.
"""
function validate(nm::Tuple{Integer,Integer})
	# CP
	
	grid  = getgrid(nm)
	A_CP  = da_LNK_cp_sz(grid, nt)
	A1_CP = copy(A_CP)
	C0_CP = CSC_LNKs_s!(A1_CP)
	C1_CP = copy(C0_CP)
	C2_CP = da_csc_LNK_cp!(C1_CP, grid; offset=1)
	
	# depth=1
	
	grid, nnts, s, nr, cfp = preparatory_multi_ps(nm, nt, 1)
	A_m1  = da_RLNK_oc_ps_sz(grid, nnts, s, cfp, nt, 1)
	C0_m1 = CSC_RLNK_si_oc_ps(A_m1, nr, s, nt, 1)
	C1_m1 = copy(C0_m1)
	C2_m1 = da_csc_RLNK_oc_ps_sz!(C1_m1, grid, nnts, s, cfp, nr, nt, 1; offset=1)
	
	# depth=2
	
	grid, nnts, s, nr, cfp = preparatory_multi_ps(nm, nt, 2)
	A_m2  = da_RLNK_oc_ps_sz(grid, nnts, s, cfp, nt, 2)
	C0_m2 = CSC_RLNK_si_oc_ps(A_m2, nr, s, nt, 2)
	C1_m2 = copy(C0_m2)
	C2_m2 = da_csc_RLNK_oc_ps_sz!(C1_m1, grid, nnts, s, cfp, nr, nt, 2; offset=1)	
		
	@warn (C0_CP-C0_m1).nzval
	@warn (C0_CP-C0_m2).nzval
	
	@warn (C2_CP-C2_m1).nzval
	@warn (C2_CP-C2_m2).nzval
	
end


#=
function validate(flag::Integer)
	
	# CP
	
	grid  = grid_[flag] #preparatory_CP(nm, nt)
	A_CP  = da_LNK_cp_sz(grid, nt)
	A1_CP = copy(A_CP)
	C0_CP = CSC_LNKs_s!(A1_CP)
	C1_CP = copy(C0_CP)
	C2_CP = da_csc_LNK_cp!(C1_CP, grid; offset=1)
	
	# depth=1
	
	grid, nnts, s, cfp, nr = gridm1_[flag], nntsm1_[flag], sm1_[flag], cfpm1_[flag], nrm1_[flag]
	A_m1  = da_RLNK_oc_ps_sz(grid, nnts, s, cfp, nt, 1)
	C0_m1 = CSC_RLNK_si_oc_ps(A_m1, nr, s, nt, 1)
	C1_m1 = copy(C0_m1)
	C2_m1 = da_csc_RLNK_oc_ps_sz!(C1_m1, grid, nnts, s, cfp, nr, nt, 1; offset=1)
	
	# depth=2
	
	grid, nnts, s, cfp, nr = gridm2_[flag], nntsm2_[flag], sm2_[flag], cfpm2_[flag], nrm2_[flag]
	A_m2  = da_RLNK_oc_ps_sz(grid, nnts, s, cfp, nt, 2)
	C0_m2 = CSC_RLNK_si_oc_ps(A_m1, nr, s, nt, 2)
	C1_m2 = copy(C0_m2)
	C2_m2 = da_csc_RLNK_oc_ps_sz!(C1_m1, grid, nnts, s, cfp, nr, nt, 2; offset=1)	
		
	@warn (C0_CP-CO_m1).nzval
	@warn (C0_CP-CO_m2).nzval
	
	@warn (C2_CP-C2_m1).nzval
	@warn (C2_CP-C2_m2).nzval
	
end
=#
