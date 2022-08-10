# VoronoiGraph
[![DOI](https://zenodo.org/badge/417525067.svg)](https://zenodo.org/badge/latestdoi/417525067)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://axsk.github.io/VoronoiGraph.jl/dev)
[![Build Status](https://github.com/axsk/VoronoiGraph.jl/workflows/CI/badge.svg)](https://github.com/axsk/VoronoiGraph.jl/actions)
[![codecov](https://codecov.io/gh/axsk/VoronoiGraph.jl/branch/main/graph/badge.svg?token=OYHZKYOE2H)](https://codecov.io/gh/axsk/VoronoiGraph.jl)


This Package implements a variation of the Voronoi Graph Traversal algorithm by Polianskii and Pokorny [\[1\]](https://dl.acm.org/doi/10.1145/3394486.3403266).
It constructs a [Voronoi Diagram](https://en.wikipedia.org/wiki/Voronoi_diagram) from a set of points by performing a random walk on the graph of the vertices of the diagram.
Unlike many other Voronoi implementations this algorithm is not limited to 2 or 3 dimensions and promises good performance even in higher dimensions.

## Usage
We can compute the Voronoi diagram with a simple call of `voronoi`
```julia
julia> data = rand(4, 100)  # 100 points in 4D space
julia> v, P = voronoi(data)
```
which returns the vertices `v::Dict`. `keys(v)` returns the simplicial complex of the diagram,
wheras `v[xs]` returns the coordinates of the vertex inbetween the generators `data[xs]`.
Additionally `P` contains the data in a vector-of-vectors format, used for further computations.

It also exports the random walk variant (returning only a subset of vertices):
```julia
julia> v, P = voronoi_random(data, 1000)  # perform 1000 iterations of the random walk
```

## Area / Volume computation
```julia
julia> A,V = volumes(v, P)
```
computes the (deterministic) areas of the boundaries of neighbouring cells (as sparse array `A`)
as well as the volume of the cells themselves (vector `V`) by falling back onto the Polyhedra.jl volume computation.


## Monte Carlo
Combining the raycasting approach with Monte Carlo estimates we can approximate the areas and volumes effectively:
```julia
julia> A, V = mc_volumes(P, 1000)  # cast 1000 Monte Carlo rays per cell
```

If the simplicial complex of vertices is already known we can speed up the process:
```julia
julia> A, V = mc_volumes(v, P, 1000)  # use the neighbourhood infromation contained in v
```

We furthermore can integrate any function `f` over a cell `i` and its boundaries:
```julia
julia> y, δy, V, A = mc_integrate(x->x^2, 1, P, 100, 10) # integrate cell 1 with 100 boundary and 100*10 volume samples
```
Here `y` and the vector `δy` contain the integrals over the cell and its boundaries. V and A get computed as a byproduct.

## References
[\[1\]](https://dl.acm.org/doi/10.1145/3394486.3403266) V. Polianskii, F. T. Pokorny - Voronoi Graph Traversal in High Dimensions with Applications to Topological Data Analysis and Piecewise Linear Interpolation (2020, Proceedings of the 26th ACM SIGKDD International Conference on Knowledge Discovery & Data Mining)
