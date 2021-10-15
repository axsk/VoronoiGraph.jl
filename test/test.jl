

function test_volume(npoints = 100, dim = 3)
	x = rand(dim, npoints)
	@time v, p = voronoi_random(x, 10000)
	@time a = adjacency(v)
	@time av = [boundary_area_verts(a, v, p) for a in a]
	@time ae = [boundary_area_edges(a, p) for a in a]
	delta = filter(isfinite, ae - av)
	quota = sum(delta .< 1e-6) / length(delta)
	@assert quota > 2/3
	@assert length(delta) / length(av) > 2/3
end


function test_random(d=6, n=200)
	x = rand(d, n)
	@time v,p = voronoi_random(rand(d,n), 10_000_000; maxstuck=100_000);
	@show length(v)
	a = adjacency(v)
	println("avg. no. of neighbors: ", length(a)/n)
	println("avg. no. of vertices per face: ", mean(length.(values(a))))
	@time c = connectivity_matrix(v, p)
	v,p,c
end


function test_grid(n=5, iter=10000)
	plot(legend=false);
	x = hcat(hexgrid(n)...)
	x .+= randn(2,n*n) .* 0.01
	@time v, P  = voronoi_random(x, iter)

	v = Dict(filter(collect(v)) do (k,v)
		norm(v) < 10
		end)

	#c = extractconn(v)
	@time A, Vs = connectivity_matrix(v, P)

	AA = map(x->x>.0, A)
	plot_connectivity!(AA .* 2, P)
	scatter!(eachrow(hcat(values(v)...))...)
	xlims!(1,6); ylims!(0,5)
end

#=
function benchmark(n=100, d=6, iter=100, particles=10)
	x = rand(d, n)
	@benchmark voronoi($x, $iter, $particles)
end
=#


function hexgrid(n)
	P = []
	for i in 1:n
		for j in 1:n
			x = i + j/2
			y = j * sqrt(3) / 2
			push!(P, [x,y])
		end
	end
	P
end

function tests()
	test_volume()
	test_random()
	test_grid()
end
