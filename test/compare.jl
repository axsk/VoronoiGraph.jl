using VoronoiCells
using VoronoiCells.GeometryBasics
using VoronoiGraph
using BenchmarkTools
using Test

function compare_area()
    data = rand(2,1000)

    rect = Rectangle(Point2(-1000., -1000.), Point2(1000., 1000))
    tess = voronoicells(data[1,:], data[2,:], rect)
    area = voronoiarea(tess) |> sort

    v, P = VoronoiGraph.voronoi(data)
    _, A = VoronoiGraph.area_volume(v, P)

    A = sort(A)

    small = A .< 1

    println("comparing areas on $(sum(small)) out of $(length(small)) cells")
    #@show maximum(abs.(A[small] - area[small]))
    @test ≈(A[small], area[small], rtol=1e-8)
end

function benchmark(n=1000)
    data = rand(2,n)

    rect = Rectangle(Point2(-1000., -1000.), Point2(1000., 1000))
    tess = voronoicells(data[1,:], data[2,:], rect)
    @show @benchmark VoronoiCells.voronoicells($data[1,:], $data[2,:], $rect)
    @show @benchmark VoronoiCells.voronoiarea($tess)

    v, P = VoronoiGraph.voronoi(data)
    @show @benchmark VoronoiGraph.voronoi($data)
    @show @benchmark VoronoiGraph.area_volume($v, $P)
    return nothing
end
