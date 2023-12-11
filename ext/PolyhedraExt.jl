module PolyhedraExt

using Polyhedra: HalfSpace, polyhedron, volume, hasrays, hrep, vrep
using VoronoiGraph: Sigma, Vertices, transformation
using LinearAlgebra
import VoronoiGraph: boundary_area

""" given two generators `g1`, `g2`, `vertices` of their common boundary
(in generator representation) and the list of all generators, compute the boundary volume.
It works by constructing a polytope from halfspaces between the generators of the boundary
vertices. We use an affine transformation to get rid of the dimension in direction g1->g2
"""
function boundary_area(g1::Int, g2::Int, vertices::AbstractVector{<:Sigma}, generators)

    length(generators[1]) == 1 && return 1

    A = generators[g1]
    B = generators[g2]
    transform = transformation(A, B)
    A = transform(A)
    gen_inds = unique!(sort!(reduce(vcat, vertices)))

    halfspaces = []
    for gen_ind in gen_inds
        gen_ind in (g1, g2) && continue
        B = transform(generators[gen_ind])
        u = normalize(B - A)
        b = dot(u, (A + B) / 2)
        hf = HalfSpace(u[2:end], b)
        push!(halfspaces, hf)
    end
    halfspaces = [h for h in halfspaces]

    # Performance with different libraries on random 6x100 data:
    # nolib: 30 mins
    # CDDLib: 40 mins
    # QHull: 5h (with warnings about not affine polyhedron using a solver)
    poly = polyhedron(hrep(halfspaces))
    local vol
    try
        vol = hasrays(poly) ? Inf : volume(poly)
    catch e
        if isa(e, AssertionError)
            @warn e
            vol = NaN
        else
            rethrow(e)
        end
    end

    return vol
end

""" similar to `boundary_area`, however uses a vector representation and is slower """
function boundary_area_vrep(g1::Int, g2::Int, inds::AbstractVector{<:Sigma}, vertices::Vertices, generators)
    A = generators[g1]
    B = generators[g2]
    dim = length(A)
    vertex_coords = map(i -> vertices[i], inds)
    push!(vertex_coords, A)  # Append one Voronoi center for full-dimensional volume.

    poly = polyhedron(vrep(vertex_coords), QHull.Library())
    vol = hasrays(poly) ? Inf : volume(poly)

    h = norm(B - A) / 2
    A = dim * vol / h
    return A, h, vol
end

end