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
