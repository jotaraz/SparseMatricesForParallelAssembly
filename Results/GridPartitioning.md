`function splitted_up!(grid::ExtendableGrid, nt::Int, numt::Int, conv, sep_par, build_par)`

200 x 401 grid:

splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_seq!, true, true)
226 ms | 350 MB

splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_seq!, false, true)
252 ms | 250 MB

splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_seq!, true, false)
218 ms | 223 MB

splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_seq!, false, false)
254 ms | 223 MB

splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_half!, true, true)
231 ms | 363 MB

splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_par!, true, true)
235 ms | 363 MB

200 x 401:


@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_seq!, false, false)
266 ms | 223 MB
@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_seq!, true, false)
244 ms | 223 MB

@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_seq!, false, true)
286 ms | 350 MB
@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_seq!, true, true)
259 ms | 350 MB

@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_half!, false, true)
288 ms | 364 MB
@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_half!, true, true)
262 ms | 364 MB

@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_par!, false, true)
288 ms | 364 MB
@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_then_CSC_par!, true, true)
265 ms | 364 MB

@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_directly_to_CSC_seq, false, true)
284 ms | 462 MB
@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_directly_to_CSC_seq, true, true)
266 ms | 462 MB

@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_directly_to_CSC_half, false, true)
261 ms | 444 MB
@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_directly_to_CSC_half, true, true)
232 ms | 444 MB

@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_directly_to_CSC_par, false, true)
254 ms | 461 MB
@warn @benchmark splitted_up!(grid1, 8, 8, Add_NonCSC_directly_to_CSC_par, true, true)
228 ms | 461 MB

800 x 401:



