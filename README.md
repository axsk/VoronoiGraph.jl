# VoronoiGraph

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://axsk.github.io/VoronoiGraph.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://axsk.github.io/VoronoiGraph.jl/dev)
[![Build Status](https://github.com/axsk/VoronoiGraph.jl/workflows/CI/badge.svg)](https://github.com/axsk/VoronoiGraph.jl/actions)

This Package implements a variation of the Voronoi Graph Traversal algorithm by Polianskii and Pokorny [\[1\]](https://dl.acm.org/doi/10.1145/3394486.3403266).
It constructs a [Voronoi Diagram](https://en.wikipedia.org/wiki/Voronoi_diagram) from a set of points by performing a random walk on the graph of the vertices of the diagram.
Unlike many other Voronoi implementations this algorithm is not limited to 2 or 3 dimensions and promises good performance even in higher dimensions.

## Usage
The main functionality is given by:
```julia
julia> data = rand(4, 100)  # 100 points in 4D space
julia> v, P = voronoi(data)
```
which returns the vertices (`v`) of the complete tesselation.

It also exports the random walk variant (returning only a subset of vertices):
```julia
julia> v, P = voronoi_random(data, 1000)  # perform 1000 iterations of the random walk
```

## Area / Volume computation
The function `area_volume` computes the (exact) areas of the boundaries of neighbouring cells
as well as the volume of the cells themselves by falling back onto the Polyhedra.jl volume computation.
This can be useful for finite volume discretizations, e.g. in the Sqra.jl (TBP) package.

Alternatively the function `mc_volumes(xs::Points)` computes an Monte-Carlo estimate of the
areas and volumes and scales better in higher dimensions.

## References
[\[1\]](https://dl.acm.org/doi/10.1145/3394486.3403266) V. Polianskii, F. T. Pokorny - Voronoi Graph Traversal in High Dimensions with Applications to Topological Data Analysis and Piecewise Linear Interpolation (2020, Proceedings of the 26th ACM SIGKDD International Conference on Knowledge Discovery & Data Mining)
