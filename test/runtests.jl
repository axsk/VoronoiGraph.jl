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
seed = rand(UInt)
@show seed
Random.seed!(seed)

# necessary condition for the voronoi diagram to be correct
function test_equidistance(verts, xs)
    allsame(x) = all(y -> y ≈ first(x), x)
    for (sig, v) in verts
        dists = [norm(xs[s] - v) for s in sig]
        relerr = maximum(abs, 1 .- dists ./ dists[1])
        if relerr > 1e-5
            @show dists
            return false
        end
    end
    return true
end

bigdata = [rand(2, 10000), rand(4,1000), rand(6,200), rand(8,100)]
smalldata = [rand(2, 1000), rand(3, 1000), rand(4,100)]

@testset "VoronoiGraph.jl" begin

    @testset "Equidistance" begin
        for data in bigdata
            @time verts, xs = voronoi(data)
            @test test_equidistance(verts, xs)

            vrand, xs = voronoi_random(data, 1000)
            @test test_equidistance(vrand, xs)

            @test issubset(keys(vrand), keys(verts))
        end
    end

    @testset "Raycast" begin
        for data in smalldata
            data = VoronoiGraph.vecvec(data)
            search = VoronoiGraph.RaycastCompare(data)
            @test_nowarn voronoi(data, search)
        end
    end

    @testset "Area computation" begin
        data = rand(2,1000)

        rect = Rectangle(Point2(-1000., -1000.), Point2(1000., 1000))
        tess = voronoicells(data[1,:], data[2,:], rect)
        area = voronoiarea(tess) |> sort

        v, P = VoronoiGraph.voronoi(data)
        _, A = VoronoiGraph.volumes(v, P)

        A = sort(A)

        # the larger areas lie at the boundary and are computed differently
        small = A .< 1
        #println("comparing areas on $(sum(small)) out of $(length(small)) cells")

        @test ≈(A[small], area[small], rtol=1e-8)
    end

    @testset "Monte Carlo Volumes" begin
        x = rand(2,20)
        v, xs = voronoi(x)
        A, V = volumes(v, xs)  # exact reference volumes

        Am, Vm = mc_volumes(xs, 10_000)
        err = std(filter(isfinite, (Am-A)./A))
        @test err < 0.2

        Am, Vm = mc_volumes(v, xs, 10_000)
        err = std(filter(isfinite, (Am-A)./A))
        @test err < 0.2
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
    @show @benchmark VoronoiGraph.volumes($v, $P)
    return nothing
end
