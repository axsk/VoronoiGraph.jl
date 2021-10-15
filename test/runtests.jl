using VoronoiGraph
using Test

@testset "VoronoiGraph.jl" begin
    verts, xs = voronoi_random(rand(2,100), 1000)
    verts, xs = voronoi(rand(2,100))
    area, vol = area_volume(verts, xs)
    adj = adjacency(verts)
end

@test_skip include("test.jl")

# TODO: turn this into a test
# begin # check if new point is equidistant to its generators
# 	rr = r + ts*u
# 	diffs = [sum(abs2, rr.-s) for s in tau]
#   allapprox(x) = all(isapprox(x[1], y) for y in x)
# 	!allapprox(diffs) && error()
# end
