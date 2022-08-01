module VoronoiGraph

using LinearAlgebra
using NearestNeighbors
using Polyhedra
using ProgressMeter
using SparseArrays
using StaticArrays
using RecipesBase

include("voronoi.jl")
include("raycast.jl")
include("volume.jl")
include("plot.jl")
include("montecarlo.jl")

export voronoi, voronoi_random
export volumes
export mc_volumes, mc_integrate

for i in 1:10
    precompile(voronoi, (Vector{SVector{i, Float64}},))
end

end
