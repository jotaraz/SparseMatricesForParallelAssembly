### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ b3235fe9-d52d-4896-933a-6dbb8bbb38e3
import Pkg; Pkg.add("ExtendableSparse"); Pkg.add("BenchmarkTools"); Pkg.add("Plots");

# ╔═╡ c054a39e-603f-11ee-3b4a-4bdc3f5c3309
begin
	using ExtendableSparse
	using BenchmarkTools
	using Plots
	using Base.Threads
	using SparseArrays
	using LinearAlgebra
	using SparseMatrixDicts
end

# ╔═╡ a319d376-55f5-4966-ae11-62cb6f6004b7




# ╔═╡ eb47422d-204d-4501-bea7-2bf9a76072cb
function dict(n, mpt, nt)
	As = [SparseMatrixDict(n, n) for i=1:nt]
	
	@threads :static for tid=1:nt
		for i=1:mpt
			As[tid][rand(1:n), rand(1:n)] += 1.0
		end
	end		
end

# ╔═╡ 5ce9735b-f3f6-4484-be5e-18815a046d03
function LNK(n, mpt, nt)
	As = [SparseMatrixLNK{Float64, Int64}(n, n) for i=1:nt]
	
	@threads :static for tid=1:nt
		for i=1:mpt
			As[tid][rand(1:n), rand(1:n)] += 1.0
		end
	end	
	
end

# ╔═╡ 9e578b8d-c775-4cd5-b30f-6a6a3a307c31
function compare(fcts, labels, n, ms, nt, ss)
	plot(;dpi=460)
	for (j,fct) in enumerate(fcts)
		times = 0.0*ms
		for (i,m) in enumerate(ms)
			@warn j, m
			for _=1:ss
				times[i] += @elapsed fct(n, m, nt)
 			end
		end

		plot!(ms, times/ss, label=labels[j])
		
	end
	plot!(xscale=:log, yscale=:log, legend=:topleft)
end

# ╔═╡ 233065b6-93b8-4de6-8388-985650eee409
compare([dict, LNK], ["dict", "LNK"], 10_000, [500, 1000, 2000, 5000, 7500, 10000, 13000, 17000, 20000, 40000, 50000, 100000, 200000], 2, 100)

# ╔═╡ Cell order:
# ╠═b3235fe9-d52d-4896-933a-6dbb8bbb38e3
# ╠═c054a39e-603f-11ee-3b4a-4bdc3f5c3309
# ╠═a319d376-55f5-4966-ae11-62cb6f6004b7
# ╠═eb47422d-204d-4501-bea7-2bf9a76072cb
# ╠═5ce9735b-f3f6-4484-be5e-18815a046d03
# ╠═9e578b8d-c775-4cd5-b30f-6a6a3a307c31
# ╠═233065b6-93b8-4de6-8388-985650eee409
