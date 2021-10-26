using Test
using VoronoiGraph

using BenchmarkTools
using LinearAlgebra
using Random
using RecipesBase
using Statistics
using VoronoiCells
using VoronoiCells.GeometryBasics

# nice to have if tests dont pass
@show Random.seed!()

# necessary condition for the voronoi diagram to be correct
function test_equidistance(verts, xs)
    allsame(x) = all(y -> y ≈ first(x), x)
    for (sig, v) in verts
        dists = [norm(xs[s] - v) for s in sig]
        !allsame(dists) && return false
    end
    return true
end


@testset "VoronoiGraph.jl" begin

    @testset "Equidistance" begin
        for (verts,xs) in [
            voronoi(rand(2,100)),
            voronoi(rand(3,100)),
            voronoi_random(rand(6,100), 1000)]

            @test test_equidistance(verts, xs)
        end
    end

    @testset "Area / Volume" begin
        data = rand(2,1000)

        rect = Rectangle(Point2(-1000., -1000.), Point2(1000., 1000))
        tess = voronoicells(data[1,:], data[2,:], rect)
        area = voronoiarea(tess) |> sort

        v, P = VoronoiGraph.voronoi(data)
        _, A = VoronoiGraph.area_volume(v, P)

        A = sort(A)
        small = A .< 1

        println("comparing areas on $(sum(small)) out of $(length(small)) cells")

        @test ≈(A[small], area[small], rtol=1e-8)
    end

    @testset "Monte Carlo Volumes" begin
        function std_error_mc_volumes(n=20, d=2, mc=10000)
            x = rand(d,n)
            v, xs = voronoi(x)
            A, V = area_volume(v, xs)
            Am, Vm = mc_volumes(xs, mc)
            std(filter(isfinite, (Am-A)./A))
        end

        @test std_error_mc_volumes() < 0.2
    end

    @testset "Plot recipe" begin
        verts, xs = voronoi(rand(2,20))
        p = VoronoiGraph.VoronoiPlot((verts, xs))

        # this is broken without `using Plots`, but works with iter
        # however we dont want to require Plots for testing, so accepting this as broken
        @test_broken RecipesBase.apply_recipe(Dict{Symbol, Any}(), p)
    end
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