

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
    gen_inds = unique!(sort!(reduce(vcat,vertices)))

    halfspaces = []
    for gen_ind in gen_inds
        gen_ind in (g1,g2) && continue
        B = transform(generators[gen_ind])
        u = normalize(B-A)
        b = dot(u, (A+B)/2)
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

""" affine transformation rotatinig and translating such that the boundary is aligned with
the first dimension. A->B will be mapped to const*[1,0,0,...] and (A+B)/2 to [0,0,...] """
function transformation(A, B)
    R = diagm(ones(length(A)))
    R[:,1] = B-A
    R = collect(inv(qr(R).Q))  # Shortcut for Gram-Schmidt orthogonalization.
    t = R*(A+B)/2
    transform(x) = R*x - t
    return transform
end

function transformation(A::SVector{N, T}, B::SVector{N, T}) where {N, T}
    R = MMatrix{N,N,T}(I)
    R[:,1] = B-A
    R = inv(qr(R).Q)
    t = R * (A+B)/2
    transform(x) = R*x - t
end

""" build the connectivity matrix for the SQRA from adjacency and boundary information """
function volumes(vertices, P::AbstractVector)
    dim = length(P[1])
    conns = adjacency(vertices)

    I = Int[]
    J = Int[]
    As = Float64[]
    Vs = zeros(length(P))

    #=@showprogress 1 "Voronoi volumes "=# for ((g1,g2), sigs) in conns
        a = boundary_area(g1, g2, sigs, P)
        h = norm(P[g1] - P[g2]) / 2
        v = a * h / dim  # Volume computation
        push!(I, g1)
        push!(J, g2)
        push!(As, a)
        Vs[g1] += v
        Vs[g2] += v
    end
    A = sparse(I, J, As, length(P), length(P))
    A = A + A'

    return A, Vs
end

""" given vertices in generator-coordinates,
collect the verts belonging to generator pairs, i.e. boundary vertices """
function adjacency(v::Vertices)
    conns = Dict{Tuple{Int,Int}, Vector{Vector{Int}}}()

    for (sig, r) in v
        for a in sig
            for b in sig
                a >= b && continue
                v = get!(conns, (a,b), [])
                push!(v, sig)
            end
        end
    end

    return conns
end


""" similar to `boundary_area`, however uses a vector representation and is slower """
function boundary_area_vrep(g1::Int, g2::Int, inds::AbstractVector{<:Sigma}, vertices::Vertices, generators)
    A = generators[g1]
    B = generators[g2]
    dim = length(A)
    vertex_coords = map(i->vertices[i], inds)
    push!(vertex_coords, A)  # Append one Voronoi center for full-dimensional volume.

    poly = polyhedron(vrep(vertex_coords), QHull.Library())
    vol = hasrays(poly) ? Inf : volume(poly)

    h = norm(B - A) / 2
    A = dim * vol / h
    return A, h, vol
end

## lets try computing the volumes ourselves, following martins idea of simplex decomposition.
# We subdivide each delauney simplex into the subsimplices with the central voronoi vertex
# instead of each of the voronoi generators as subsimplex-vertex. (1)
# This smaller simplex contributes equal volume to all adjacent cells. (2)
# We further can compute the areas by means of the pyramid formula

function myvolumes(vertices::Vertices, P::Points)
    dim = length(P[1])

    As = Float64[]
    Vs = zeros(length(P))

    for v in vertices
        σ, x = v
        X = reduce(hcat, P[σ]) .- x
        for i in 1:dim+1
            mask = 1:dim+1 .!= i
            A = X[:, mask]
            #if dot(X[:,i], A[:,1])  # vertex i
            V = abs(det(A)) / factorial(dim)  # (1)
            Vs[σ[mask]] .+= V / dim  # (2)
        end
    end
    return Vs
end

# this does not work: sometimes the voronoi vertex does not lie inside but outside of
# its generator's simplex. this leads to missattribution of volumes.
# seems not so easy to handle as its nonlocal (to the vertex) information
