using VoronoiGraph
using Test

@testset "VoronoiGraph.jl" begin
	verts, xs = voronoi_random(rand(2,100), 1000)
    verts, xs = voronoi(rand(2,100))
	area, vol = area_volume(verts, xs)
	adj = adjacency(verts)
end

@test_skip include("test.jl")
