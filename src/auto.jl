
#path = "/home/johannes/Nextcloud/Documents/Uni/VIII/WIAS/juliaCode(anfang)/para/somezeros/"#"/Home/.../"
#pre_name = "/home/johannes/Nextcloud/Documents/Uni/VIII/WIAS/juliaCode(anfang)/para/somezeros/LenovoT480s/"

path = "/Home/guests/taraz/somezeros/" #"/Home/.../"
pre_name = path*"escher-04/"


include(path*"benchmark.jl")

#2d
ns = [400, 800, 1200, 1600, 2000]
for n in ns
	add_a_grid((n, n+1))
end

#3d
ns = [32, 64, 96, 128]
for n in ns
	add_a_grid((n, n, n+1))
end

run_auto_nowrite(;ss=2)

run_auto_allocs()

run_auto(1, 10)
run_auto(2, 10)
run_auto(3, 10)
