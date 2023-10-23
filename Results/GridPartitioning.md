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

Set1:

seq! f f  | 1.514 s | 895.60 MB
seq! t f  | 1.403 s | 895.60 MB
seq! f t  | 1.331 s | 1.37 GB
seq! t t  | 1.369 s | 1.37 GB
half! f t | 1.504 s | 1.43 GB
half! t t | 1.236 s | 1.43 GB
par! f t  | 1.503 s | 1.43 GB
par! t t  | 1.222 s | 1.43 GB
seq f t   | 1.596 s | 1.81 GB
seq t t   | 1.251 s | 1.81 GB
half f t  | 1.187 s | 1.74 GB
half t t  | 1.025 s | 1.74 GB
par f t   | 1.147 s | 1.81 GB
par t t   | 1.037 s | 1.81 GB

-------------------------------------------------
Set2:

seq! f f  | 1.338 s | 895.60 MB
seq! t f  | 1.160 s | 895.60 MB
seq! f t  | 1.401 s | 1.37 GB
seq! t t  | 1.221 s | 1.37 GB
half! f t | 1.508 s | 1.43 GB
half! t t | 1.285 s | 1.43 GB
par! f t  | 1.427 s | 1.43 GB
par! t t  | 1.295 s | 1.43 GB
seq f t   | 1.606 s | 1.81 GB
seq t t   | 1.494 s | 1.81 GB
half f t  | 1.530 s | 1.74 GB
half t t  | 1.324 s | 1.74 GB
par f t   | 1.426 s | 1.81 GB
par t t   | 1.306 s | 1.81 GB


--------------------------------------------------------


800 x 801


Set 1

seq! f f  | 3.228 s | 1.54 GB
seq! t f  | 2.985 s | 1.75 GB
seq! f t  | 3.408 s | 2.74 GB
seq! t t  | 2.675 s | 2.74 GB
half! f t | 3.493 s | 2.85 GB
half! t t | 2.926 s | 2.85 GB
par! f t  | 3.680 s | 2.85 GB
par! t t  | 3.173 s | 2.85 GB
seq f t   | 3.937 s | 3.62 GB
seq t t   | 3.093 s | 3.62 GB
half f t  | 3.140 s | 3.48 GB
half t t  | 2.569 s | 3.48 GB
par f t   | 2.852 s | 3.61 GB
par t t   | 2.574 s | 3.61 GiB
 
 
 Set 2
 
seq! f f  | 3.428 s | 1.75 GB
seq! t f  | 3.025 s | 1.75 GB
seq! f t  | 3.619 s | 2.74 GB
seq! t t  | 3.274 s | 2.74 GB
half! f t | 3.714 s | 2.85 GB
half! t t | 3.360 s | 2.85 GB
par! f t  | 3.869 s | 2.85 GB
par! t t  | 3.445 s | 2.85 GB
seq f t   | 3.801 s | 3.62 GB
seq t t   | 3.457 s | 3.62 GB
half f t  | 3.301 s | 3.48 GB
half t t  | 3.098 s | 3.48 GB
par f t   | 3.153 s | 3.61 GB
par t t   | 2.898 s | 3.61 GB



-------------------------------------------------------------------


800 x 401:

----------------------------------------------------------

seq! f f
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.213 s …   1.385 s  ┊ GC (min … max): 4.53% … 7.83%
 Time  (median):     1.315 s              ┊ GC (median):    7.86%
 Time  (mean ± σ):   1.307 s ± 83.120 ms  ┊ GC (mean ± σ):  7.48% ± 2.07%

  █               █                                  █    █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁█ ▁
  1.21 s         Histogram: frequency by time        1.39 s <

 Memory estimate: 895.60 MiB, allocs estimate: 5113659.
seq! t f
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.053 s …   1.172 s  ┊ GC (min … max): 4.91% … 9.45%
 Time  (median):     1.108 s              ┊ GC (median):    7.82%
 Time  (mean ± σ):   1.112 s ± 42.845 ms  ┊ GC (mean ± σ):  7.61% ± 1.66%

  █                      █ █        █                     █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁█▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.05 s         Histogram: frequency by time        1.17 s <

 Memory estimate: 895.60 MiB, allocs estimate: 5113710.
seq! f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.303 s …   1.391 s  ┊ GC (min … max): 3.18% … 7.42%
 Time  (median):     1.354 s              ┊ GC (median):    6.14%
 Time  (mean ± σ):   1.351 s ± 45.177 ms  ┊ GC (mean ± σ):  5.76% ± 1.90%

  █          █                                          █ █  
  █▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁█ ▁
  1.3 s          Histogram: frequency by time        1.39 s <

 Memory estimate: 1.37 GiB, allocs estimate: 5113858.
seq! t t
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.195 s …   1.347 s  ┊ GC (min … max): 3.91% … 10.93%
 Time  (median):     1.242 s              ┊ GC (median):    8.29%
 Time  (mean ± σ):   1.261 s ± 70.743 ms  ┊ GC (mean ± σ):  7.78% ±  2.56%

  █                ▁                              ▁       ▁  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁█ ▁
  1.2 s          Histogram: frequency by time        1.35 s <

 Memory estimate: 1.37 GiB, allocs estimate: 5113909.
half! f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.368 s …   1.487 s  ┊ GC (min … max): 3.25% … 7.52%
 Time  (median):     1.474 s              ┊ GC (median):    7.63%
 Time  (mean ± σ):   1.451 s ± 55.889 ms  ┊ GC (mean ± σ):  6.71% ± 2.28%

  █                                             █       █ █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁█▁█ ▁
  1.37 s         Histogram: frequency by time        1.49 s <

 Memory estimate: 1.43 GiB, allocs estimate: 5113911.
half! t t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.285 s …   1.309 s  ┊ GC (min … max): 3.77% … 8.90%
 Time  (median):     1.294 s              ┊ GC (median):    8.36%
 Time  (mean ± σ):   1.295 s ± 11.325 ms  ┊ GC (mean ± σ):  7.36% ± 2.40%

  █   █                                █                  █  
  █▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.28 s         Histogram: frequency by time        1.31 s <

 Memory estimate: 1.43 GiB, allocs estimate: 5113962.
par! f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.431 s …    1.677 s  ┊ GC (min … max): 2.82% … 6.68%
 Time  (median):     1.487 s               ┊ GC (median):    7.46%
 Time  (mean ± σ):   1.520 s ± 107.517 ms  ┊ GC (mean ± σ):  6.25% ± 2.32%

  ▁            █                                           ▁  
  █▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.43 s         Histogram: frequency by time         1.68 s <

 Memory estimate: 1.43 GiB, allocs estimate: 5113961.
par! t t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.302 s …   1.364 s  ┊ GC (min … max): 7.95% … 9.35%
 Time  (median):     1.346 s              ┊ GC (median):    8.23%
 Time  (mean ± σ):   1.339 s ± 27.372 ms  ┊ GC (mean ± σ):  7.40% ± 2.64%

  █                                █             █        █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁█ ▁
  1.3 s          Histogram: frequency by time        1.36 s <

 Memory estimate: 1.43 GiB, allocs estimate: 5114013.
seq f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.345 s …    1.980 s  ┊ GC (min … max):  4.49% … 24.61%
 Time  (median):     1.766 s               ┊ GC (median):    19.72%
 Time  (mean ± σ):   1.714 s ± 276.709 ms  ┊ GC (mean ± σ):  19.05% ± 11.40%

  █                            █                █          █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁█ ▁
  1.35 s         Histogram: frequency by time         1.98 s <

 Memory estimate: 1.81 GiB, allocs estimate: 5113971.
seq t t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.404 s …    1.702 s  ┊ GC (min … max):  4.27% … 27.58%
 Time  (median):     1.515 s               ┊ GC (median):    11.38%
 Time  (mean ± σ):   1.534 s ± 131.931 ms  ┊ GC (mean ± σ):  14.25% ± 10.22%

  █         █                     █                        █  
  █▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.4 s          Histogram: frequency by time          1.7 s <

 Memory estimate: 1.81 GiB, allocs estimate: 5114023.
half f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.217 s …   1.282 s  ┊ GC (min … max): 6.91% … 6.71%
 Time  (median):     1.258 s              ┊ GC (median):    6.76%
 Time  (mean ± σ):   1.254 s ± 29.826 ms  ┊ GC (mean ± σ):  6.49% ± 2.34%

  █                     █                           █     █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁█ ▁
  1.22 s         Histogram: frequency by time        1.28 s <

 Memory estimate: 1.74 GiB, allocs estimate: 5114100.
half t t
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.078 s …   1.200 s  ┊ GC (min … max): 9.41% … 9.14%
 Time  (median):     1.106 s              ┊ GC (median):    9.17%
 Time  (mean ± σ):   1.118 s ± 48.276 ms  ┊ GC (mean ± σ):  8.12% ± 2.10%

  █    █       █   █                                      █  
  █▁▁▁▁█▁▁▁▁▁▁▁█▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.08 s         Histogram: frequency by time         1.2 s <

 Memory estimate: 1.74 GiB, allocs estimate: 5114151.
par f t
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.182 s …   1.293 s  ┊ GC (min … max): 3.96% … 11.03%
 Time  (median):     1.218 s              ┊ GC (median):    7.15%
 Time  (mean ± σ):   1.224 s ± 44.502 ms  ┊ GC (mean ± σ):  7.35% ±  3.07%

  █   █             █          █                          █  
  █▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.18 s         Histogram: frequency by time        1.29 s <

 Memory estimate: 1.81 GiB, allocs estimate: 5114114.
par t t
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.025 s …   1.177 s  ┊ GC (min … max): 5.20% … 7.89%
 Time  (median):     1.107 s              ┊ GC (median):    8.39%
 Time  (mean ± σ):   1.090 s ± 64.398 ms  ┊ GC (mean ± σ):  8.29% ± 2.59%

  █                             ▁ ▁                       ▁  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.02 s         Histogram: frequency by time        1.18 s <

 Memory estimate: 1.81 GiB, allocs estimate: 5114164.

---------------------------------------------------

Set2:
seq! f f
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.292 s …   1.372 s  ┊ GC (min … max): 4.27% … 10.19%
 Time  (median):     1.344 s              ┊ GC (median):    7.74%
 Time  (mean ± σ):   1.338 s ± 36.283 ms  ┊ GC (mean ± σ):  7.53% ±  2.45%

  █                       █                        █      █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁█ ▁
  1.29 s         Histogram: frequency by time        1.37 s <

 Memory estimate: 895.60 MiB, allocs estimate: 5113659.
seq! t f
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.104 s …   1.246 s  ┊ GC (min … max): 4.16% … 9.58%
 Time  (median):     1.130 s              ┊ GC (median):    7.82%
 Time  (mean ± σ):   1.160 s ± 57.665 ms  ┊ GC (mean ± σ):  7.39% ± 1.98%

  █        ██                       █                     █  
  █▁▁▁▁▁▁▁▁██▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.1 s          Histogram: frequency by time        1.25 s <

 Memory estimate: 895.60 MiB, allocs estimate: 5113710.
seq! f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.357 s …   1.463 s  ┊ GC (min … max): 3.63% … 5.75%
 Time  (median):     1.392 s              ┊ GC (median):    6.50%
 Time  (mean ± σ):   1.401 s ± 44.926 ms  ┊ GC (mean ± σ):  5.89% ± 1.62%

  █              █      █                                 █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.36 s         Histogram: frequency by time        1.46 s <

 Memory estimate: 1.37 GiB, allocs estimate: 5113858.
seq! t t
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.185 s …   1.257 s  ┊ GC (min … max): 4.52% … 8.14%
 Time  (median):     1.206 s              ┊ GC (median):    8.48%
 Time  (mean ± σ):   1.221 s ± 32.377 ms  ┊ GC (mean ± σ):  8.18% ± 2.94%

  █            █  █                                     █ █  
  █▁▁▁▁▁▁▁▁▁▁▁▁█▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁█ ▁
  1.19 s         Histogram: frequency by time        1.26 s <

 Memory estimate: 1.37 GiB, allocs estimate: 5113908.
half! f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.390 s …   1.593 s  ┊ GC (min … max): 3.18% … 7.74%
 Time  (median):     1.524 s              ┊ GC (median):    7.70%
 Time  (mean ± σ):   1.508 s ± 85.451 ms  ┊ GC (mean ± σ):  6.76% ± 2.39%

  █                                 █     █               █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.39 s         Histogram: frequency by time        1.59 s <

 Memory estimate: 1.43 GiB, allocs estimate: 5113912.
half! t t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.265 s …   1.314 s  ┊ GC (min … max): 3.78% … 9.26%
 Time  (median):     1.281 s              ┊ GC (median):    8.15%
 Time  (mean ± σ):   1.285 s ± 21.582 ms  ┊ GC (mean ± σ):  7.36% ± 2.54%

  █        █                 █                            █  
  █▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.26 s         Histogram: frequency by time        1.31 s <

 Memory estimate: 1.43 GiB, allocs estimate: 5113961.
par! f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.295 s …   1.513 s  ┊ GC (min … max): 3.18% … 8.27%
 Time  (median):     1.451 s              ┊ GC (median):    7.76%
 Time  (mean ± σ):   1.427 s ± 93.170 ms  ┊ GC (mean ± σ):  6.99% ± 2.62%

  ▁                                       █               ▁  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.29 s         Histogram: frequency by time        1.51 s <

 Memory estimate: 1.43 GiB, allocs estimate: 5113960.
par! t t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.236 s …   1.341 s  ┊ GC (min … max): 7.28% … 3.88%
 Time  (median):     1.301 s              ┊ GC (median):    7.39%
 Time  (mean ± σ):   1.295 s ± 46.031 ms  ┊ GC (mean ± σ):  7.00% ± 2.26%

  █                       █                    █          █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁█ ▁
  1.24 s         Histogram: frequency by time        1.34 s <

 Memory estimate: 1.43 GiB, allocs estimate: 5114013.
seq f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.296 s …    1.864 s  ┊ GC (min … max):  4.32% … 28.52%
 Time  (median):     1.633 s               ┊ GC (median):    21.92%
 Time  (mean ± σ):   1.606 s ± 278.725 ms  ┊ GC (mean ± σ):  21.69% ± 13.82%

  █              █                                     █   █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁█ ▁
  1.3 s          Histogram: frequency by time         1.86 s <

 Memory estimate: 1.81 GiB, allocs estimate: 5113971.
seq t t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.361 s …    1.712 s  ┊ GC (min … max):  4.75% … 32.73%
 Time  (median):     1.451 s               ┊ GC (median):    12.50%
 Time  (mean ± σ):   1.494 s ± 155.866 ms  ┊ GC (mean ± σ):  16.53% ± 12.30%

  █      █              █                                  █  
  █▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.36 s         Histogram: frequency by time         1.71 s <

 Memory estimate: 1.81 GiB, allocs estimate: 5114024.
half f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.420 s …    1.679 s  ┊ GC (min … max): 4.32% … 8.85%
 Time  (median):     1.510 s               ┊ GC (median):    6.83%
 Time  (mean ± σ):   1.530 s ± 108.789 ms  ┊ GC (mean ± σ):  6.80% ± 1.87%

  █               █      █                                 █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.42 s         Histogram: frequency by time         1.68 s <

 Memory estimate: 1.74 GiB, allocs estimate: 5114101.
half t t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.280 s …   1.356 s  ┊ GC (min … max): 4.90% … 9.57%
 Time  (median):     1.329 s              ┊ GC (median):    9.46%
 Time  (mean ± σ):   1.324 s ± 36.448 ms  ┊ GC (mean ± σ):  8.53% ± 2.47%

  █                   █                                █  █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁█ ▁
  1.28 s         Histogram: frequency by time        1.36 s <

 Memory estimate: 1.74 GiB, allocs estimate: 5114152.
par f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.385 s …   1.465 s  ┊ GC (min … max): 9.01% … 8.40%
 Time  (median):     1.428 s              ┊ GC (median):    7.74%
 Time  (mean ± σ):   1.426 s ± 42.441 ms  ┊ GC (mean ± σ):  7.37% ± 1.66%

  █     █                                               █ █  
  █▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁█ ▁
  1.38 s         Histogram: frequency by time        1.46 s <

 Memory estimate: 1.81 GiB, allocs estimate: 5114114.
par t t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.263 s …   1.385 s  ┊ GC (min … max): 10.19% … 12.34%
 Time  (median):     1.287 s              ┊ GC (median):     8.28%
 Time  (mean ± σ):   1.306 s ± 54.327 ms  ┊ GC (mean ± σ):   8.58% ±  3.36%

  █      █       █                                        █  
  █▁▁▁▁▁▁█▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.26 s         Histogram: frequency by time        1.38 s <

 Memory estimate: 1.81 GiB, allocs estimate: 5114163.

----------------------------------------------------

 Set1:

seq! f f
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.469 s …   1.577 s  ┊ GC (min … max): 8.26% … 2.97%
 Time  (median):     1.505 s              ┊ GC (median):    7.35%
 Time  (mean ± σ):   1.514 s ± 45.338 ms  ┊ GC (mean ± σ):  6.68% ± 2.70%

  █                █ █                                    █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.47 s         Histogram: frequency by time        1.58 s <

 Memory estimate: 895.60 MiB, allocs estimate: 5113659.
seq! t f
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.267 s …    1.543 s  ┊ GC (min … max): 5.20% … 6.54%
 Time  (median):     1.401 s               ┊ GC (median):    8.01%
 Time  (mean ± σ):   1.403 s ± 118.915 ms  ┊ GC (mean ± σ):  7.53% ± 1.96%

  █                 █                  █                   █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.27 s         Histogram: frequency by time         1.54 s <

 Memory estimate: 895.60 MiB, allocs estimate: 5113709.
seq! f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.270 s …   1.413 s  ┊ GC (min … max): 7.34% … 5.22%
 Time  (median):     1.320 s              ┊ GC (median):    6.33%
 Time  (mean ± σ):   1.331 s ± 71.110 ms  ┊ GC (mean ± σ):  5.78% ± 1.86%

  █                                      ▁                ▁  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.27 s         Histogram: frequency by time        1.41 s <

 Memory estimate: 1.37 GiB, allocs estimate: 5113857.
seq! t t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.358 s …   1.379 s  ┊ GC (min … max): 6.54% … 9.15%
 Time  (median):     1.370 s              ┊ GC (median):    7.58%
 Time  (mean ± σ):   1.369 s ± 10.262 ms  ┊ GC (mean ± σ):  7.03% ± 2.49%

  █             █                                    █    █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁█ ▁
  1.36 s         Histogram: frequency by time        1.38 s <

 Memory estimate: 1.37 GiB, allocs estimate: 5113909.
half! f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.424 s …   1.595 s  ┊ GC (min … max): 2.81% … 7.26%
 Time  (median):     1.499 s              ┊ GC (median):    7.66%
 Time  (mean ± σ):   1.504 s ± 75.184 ms  ┊ GC (mean ± σ):  6.98% ± 2.95%

  █             █                    █                    █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.42 s         Histogram: frequency by time         1.6 s <

 Memory estimate: 1.43 GiB, allocs estimate: 5113909.
half! t t
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.122 s …   1.293 s  ┊ GC (min … max): 4.55% … 8.50%
 Time  (median):     1.260 s              ┊ GC (median):    8.72%
 Time  (mean ± σ):   1.236 s ± 66.483 ms  ┊ GC (mean ± σ):  8.18% ± 2.52%

  █                                     █       ██        █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁██▁▁▁▁▁▁▁▁█ ▁
  1.12 s         Histogram: frequency by time        1.29 s <

 Memory estimate: 1.43 GiB, allocs estimate: 5113963.
par! f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.351 s …    1.631 s  ┊ GC (min … max): 3.29% … 7.08%
 Time  (median):     1.514 s               ┊ GC (median):    7.58%
 Time  (mean ± σ):   1.503 s ± 121.329 ms  ┊ GC (mean ± σ):  6.78% ± 2.44%

  █                       █                  █             █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.35 s         Histogram: frequency by time         1.63 s <

 Memory estimate: 1.43 GiB, allocs estimate: 5113962.
par! t t
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.160 s …   1.278 s  ┊ GC (min … max): 4.38% … 11.79%
 Time  (median):     1.233 s              ┊ GC (median):    8.47%
 Time  (mean ± σ):   1.222 s ± 43.618 ms  ┊ GC (mean ± σ):  8.28% ±  2.70%

  █                    █             ██                   █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁██▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.16 s         Histogram: frequency by time        1.28 s <

 Memory estimate: 1.43 GiB, allocs estimate: 5114014.
seq f t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.277 s …    1.840 s  ┊ GC (min … max):  3.90% … 28.20%
 Time  (median):     1.633 s               ┊ GC (median):    21.77%
 Time  (mean ± σ):   1.596 s ± 274.144 ms  ┊ GC (mean ± σ):  21.39% ± 13.58%

  █                 █                                   █  █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁█ ▁
  1.28 s         Histogram: frequency by time         1.84 s <

 Memory estimate: 1.81 GiB, allocs estimate: 5113972.
seq t t
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.146 s …   1.296 s  ┊ GC (min … max): 4.71% … 9.79%
 Time  (median):     1.281 s              ┊ GC (median):    7.75%
 Time  (mean ± σ):   1.251 s ± 71.097 ms  ┊ GC (mean ± σ):  7.58% ± 2.10%

  ▁                                            ▁          █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁█ ▁
  1.15 s         Histogram: frequency by time         1.3 s <

 Memory estimate: 1.81 GiB, allocs estimate: 5114020.
half f t
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.117 s …   1.317 s  ┊ GC (min … max): 3.76% … 9.11%
 Time  (median):     1.149 s              ┊ GC (median):    7.04%
 Time  (mean ± σ):   1.187 s ± 79.061 ms  ┊ GC (mean ± σ):  6.98% ± 1.99%

  ▁        █             ▁                                ▁  
  █▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.12 s         Histogram: frequency by time        1.32 s <

 Memory estimate: 1.74 GiB, allocs estimate: 5114099.
half t t
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  990.433 ms …   1.054 s  ┊ GC (min … max): 10.18% … 9.19%
 Time  (median):        1.029 s              ┊ GC (median):     9.72%
 Time  (mean ± σ):      1.025 s ± 23.435 ms  ┊ GC (mean ± σ):   8.77% ± 2.38%

  █                       █            █   █                 █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  990 ms          Histogram: frequency by time          1.05 s <

 Memory estimate: 1.74 GiB, allocs estimate: 5114151.
par f t
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  1.112 s …   1.177 s  ┊ GC (min … max): 4.83% … 9.22%
 Time  (median):     1.152 s              ┊ GC (median):    9.18%
 Time  (mean ± σ):   1.147 s ± 23.476 ms  ┊ GC (mean ± σ):  7.78% ± 2.21%

  █                         █       ██                    █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁██▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.11 s         Histogram: frequency by time        1.18 s <

 Memory estimate: 1.81 GiB, allocs estimate: 5114112.
par t t
BenchmarkTools.Trial: 5 samples with 1 evaluation.
 Range (min … max):  963.637 ms …   1.083 s  ┊ GC (min … max): 4.75% … 12.62%
 Time  (median):        1.046 s              ┊ GC (median):    9.39%
 Time  (mean ± σ):      1.037 s ± 50.126 ms  ┊ GC (mean ± σ):  9.04% ±  2.88%

  █                       █                █                ██  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁██ ▁
  964 ms          Histogram: frequency by time          1.08 s <

 Memory estimate: 1.81 GiB, allocs estimate: 5114164.

