var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = VoronoiGraph","category":"page"},{"location":"#VoronoiGraph","page":"Home","title":"VoronoiGraph","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for VoronoiGraph.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [VoronoiGraph]","category":"page"},{"location":"#VoronoiGraph.adjacency-Tuple{Dict{var\"#s3\", var\"#s2\"} where {var\"#s3\"<:(AbstractVector{var\"#s2\"} where var\"#s2\"<:Integer), var\"#s2\"<:(AbstractVector{T} where T<:Real)}}","page":"Home","title":"VoronoiGraph.adjacency","text":"given vertices in generator-coordinates, collect the verts belonging to generator pairs, i.e. boundary vertices \n\n\n\n\n\n","category":"method"},{"location":"#VoronoiGraph.area_volume-Tuple{Any, AbstractVector{T} where T}","page":"Home","title":"VoronoiGraph.area_volume","text":"build the connectivity matrix for the SQRA from adjacency and boundary information \n\n\n\n\n\n","category":"method"},{"location":"#VoronoiGraph.boundary_area-Tuple{Int64, Int64, AbstractVector{var\"#s5\"} where var\"#s5\"<:(AbstractVector{var\"#s2\"} where var\"#s2\"<:Integer), Any}","page":"Home","title":"VoronoiGraph.boundary_area","text":"given two generators g1, g2, vertices of their common boundary (in generator representation) and the list of all generators, compute the boundary volume. It works by constructing a polytope from halfspaces between the generators of the boundary vertices. We use an affine transformation to get rid of the dimension in direction g1->g2\n\n\n\n\n\n","category":"method"},{"location":"#VoronoiGraph.boundary_area_vrep-Tuple{Int64, Int64, AbstractVector{var\"#s6\"} where var\"#s6\"<:(AbstractVector{var\"#s2\"} where var\"#s2\"<:Integer), Dict{var\"#s3\", var\"#s2\"} where {var\"#s3\"<:(AbstractVector{var\"#s2\"} where var\"#s2\"<:Integer), var\"#s2\"<:(AbstractVector{T} where T<:Real)}, Any}","page":"Home","title":"VoronoiGraph.boundary_area_vrep","text":"similar to boundary_area, however uses a vector representation and is slower \n\n\n\n\n\n","category":"method"},{"location":"#VoronoiGraph.descent","page":"Home","title":"VoronoiGraph.descent","text":"starting at given points, run the ray shooting descent to find vertices \n\n\n\n\n\n","category":"function"},{"location":"#VoronoiGraph.explore-Tuple{Any, Any, AbstractVector{var\"#s2\"} where var\"#s2\"<:(AbstractVector{T} where T<:Real), Any}","page":"Home","title":"VoronoiGraph.explore","text":"BFS of vertices starting from S0 \n\n\n\n\n\n","category":"method"},{"location":"#VoronoiGraph.randray-Tuple{AbstractVector{var\"#s2\"} where var\"#s2\"<:(AbstractVector{T} where T<:Real)}","page":"Home","title":"VoronoiGraph.randray","text":"generate a random ray orthogonal to the subspace spanned by the given points \n\n\n\n\n\n","category":"method"},{"location":"#VoronoiGraph.raycast-Tuple{AbstractVector{var\"#s2\"} where var\"#s2\"<:Integer, AbstractVector{T} where T<:Real, AbstractVector{T} where T<:Real, AbstractVector{var\"#s2\"} where var\"#s2\"<:(AbstractVector{T} where T<:Real), VoronoiGraph.SearchBisection}","page":"Home","title":"VoronoiGraph.raycast","text":"shooting a ray in the given direction, find the next connecting point. This variant (by Poliaski, Pokorny) uses a binary search \n\n\n\n\n\n","category":"method"},{"location":"#VoronoiGraph.raycast-Tuple{AbstractVector{var\"#s2\"} where var\"#s2\"<:Integer, AbstractVector{T} where T<:Real, AbstractVector{T} where T<:Real, AbstractVector{var\"#s2\"} where var\"#s2\"<:(AbstractVector{T} where T<:Real), VoronoiGraph.SearchIncircle}","page":"Home","title":"VoronoiGraph.raycast","text":"Shooting a ray in the given direction, find the next connecting point. This variant uses an iterative NN search \n\n\n\n\n\n","category":"method"},{"location":"#VoronoiGraph.raycast-Tuple{AbstractVector{var\"#s2\"} where var\"#s2\"<:Integer, Any, Any, Any, VoronoiGraph.SearchBruteforce}","page":"Home","title":"VoronoiGraph.raycast","text":"shooting a ray in the given direction, find the next connecting point. This is the bruteforce variant, using a linear search to find the closest point \n\n\n\n\n\n","category":"method"},{"location":"#VoronoiGraph.transformation-Tuple{Any, Any}","page":"Home","title":"VoronoiGraph.transformation","text":"affine transformation rotatinig and translating such that the boundary is aligned with the first dimension. A->B will be mapped to const*[1,0,0,...] and (A+B)/2 to [0,0,...] \n\n\n\n\n\n","category":"method"},{"location":"#VoronoiGraph.voronoi-Tuple{AbstractVector{var\"#s2\"} where var\"#s2\"<:(AbstractVector{T} where T<:Real)}","page":"Home","title":"VoronoiGraph.voronoi","text":"construct the voronoi diagram from x through breadth-first search \n\n\n\n\n\n","category":"method"},{"location":"#VoronoiGraph.voronoi_random","page":"Home","title":"VoronoiGraph.voronoi_random","text":"construct a (partial) voronoi diagram from x through a random walk \n\n\n\n\n\n","category":"function"},{"location":"#VoronoiGraph.walk-Tuple{AbstractVector{var\"#s2\"} where var\"#s2\"<:Integer, AbstractVector{T} where T<:Real, Int64, AbstractVector{var\"#s2\"} where var\"#s2\"<:(AbstractVector{T} where T<:Real), Any, Int64}","page":"Home","title":"VoronoiGraph.walk","text":"starting at vertices, walk nsteps along the voronoi graph to find new vertices \n\n\n\n\n\n","category":"method"},{"location":"#VoronoiGraph.walkray-NTuple{4, Any}","page":"Home","title":"VoronoiGraph.walkray","text":"starting at vertex (v,r), return a random adjacent vertex \n\n\n\n\n\n","category":"method"},{"location":"#VoronoiGraph.walkray-Tuple{AbstractVector{var\"#s2\"} where var\"#s2\"<:Integer, AbstractVector{T} where T<:Real, AbstractVector{var\"#s2\"} where var\"#s2\"<:(AbstractVector{T} where T<:Real), Any, Any}","page":"Home","title":"VoronoiGraph.walkray","text":"find the vertex connected to v by moving away from its i-th generator \n\n\n\n\n\n","category":"method"}]
}
