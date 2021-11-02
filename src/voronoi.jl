const Point{T} = AbstractVector{T} where T<:Real
const Points = AbstractVector{<:Point}
const Sigma = AbstractVector{<:Integer}  # Encoded the vertex by the ids of its generators.
const Vertex = Tuple{<:Sigma, <: Point}
const Vertices = Dict{<:Sigma, <:Point}


voronoi(x; kwargs...) = voronoi(vecvec(x); kwargs...)

""" construct the voronoi diagram from `x` through breadth-first search """
function voronoi(xs::Points, searcher = Raycast(xs))
    sig, r = descent(xs, searcher)
    verts = explore(sig, r, xs, searcher)
    return verts::Vertices, xs
end

voronoi_random(x, args...; kwargs...) = voronoi_random(vecvec(x), args...; kwargs...)

""" construct a (partial) voronoi diagram from `x` through a random walk """
function voronoi_random(xs::Points, iter=1000; tmax=1000, maxstuck=typemax(Int))
    searcher = RaycastIncircle(KDTree(xs), tmax)
    sig, r = descent(xs, searcher)
    verts = walk(sig, r, iter, xs, searcher, maxstuck)
    return verts::Vertices, xs
end

vecvec(x::Matrix) = map(SVector{size(x,1)}, eachcol(x))
vecvec(x::Vector{<:SVector}) = x

""" starting at given points, run the ray shooting descent to find vertices """
function descent(xs::Points, searcher, start = 1) :: Vertex
    sig = [start]
    r = xs[start]
    d = length(r)

    for k in d:-1:1  # find an additional generator for each dimension
        u = randray(xs[sig])
        (tau, t) = raycast(sig, r, u, xs, searcher)
        if t == Inf
            u = -u
            (tau, t) = raycast(sig, r, u, xs, searcher)
        end
        if t == Inf
            error("Could not find a vertex in both directions of current point." *
                "Consider increasing search range (tmax)")
        end
        sig = tau
        r = r + t*u
    end

    return (sig, r)
end


""" starting at vertices, walk nsteps along the voronoi graph to find new vertices """
function walk(sig::Sigma, r::Point, nsteps::Int, xs::Points, searcher, maxstuck::Int) :: Vertices
    verts = Dict(sig => r)
    nonew = 0
    prog = Progress(maxstuck, 1, "Voronoi walk")
    progmax = 0

    for s in 1:nsteps
        nonew += 1
        progmax = max(progmax, nonew)
        ProgressMeter.update!(prog, progmax)
        sig, r = walkray(sig, r, xs, searcher)
        get!(verts, sig) do
            nonew = 0
            return r
        end
        nonew > maxstuck && break
    end

    return verts
end

""" starting at vertex (v,r), return a random adjacent vertex """
walkray(sig, r, xs, searcher) = walkray(sig, r, xs, searcher, rand(1:length(sig)))

""" find the vertex connected to `v` by moving away from its `i`-th generator """
function walkray(sig::Sigma, r::Point, xs::Points, searcher, i)
    sig_del = deleteat(sig, i)
    u = randray(xs[sig_del])
    if (u' * (xs[sig[i]] - xs[sig_del[1]])) > 0
        u = -u
    end
    sig_new, t = raycast(sig_del, r, u, xs, searcher)
    if t < Inf
        sig = sig_new
        r = r + t*u
    end
    return sig, r
end

""" BFS of vertices starting from `S0` """
function explore(sig, r, xs::Points, searcher) :: Vertices
    verts = Dict(sig=>r)
    queue = copy(verts)
    edgecount = Dict{Vector{Int64}, Int}()

    cache = true  # Caches if both points along an edge are known. Trades memory for runtime.
    hit = 0
    miss = 0
    new = 0
    unbounded = 0

    while length(queue) > 0
        (sig,r) = pop!(queue)
        for i in 1:length(sig)

            if cache && get(edgecount, deleteat(sig, i), 0) == 2
                hit += 1
                continue
            end

            sig_new, r_new = walkray(sig, r, xs, searcher, i)

            if sig_new == sig
                unbounded += 1
                continue
            end

            if !haskey(verts, sig_new)
                push!(queue, sig_new => r_new)
                push!(verts, sig_new => r_new)
                if cache
                    for j in 1:length(sig_new)
                        edge = deleteat(sig_new, j)
                        edgecount[edge] = get(edgecount, edge, 0) + 1
                    end
                end
                new += 1
            else
                miss += 1
            end
        end

    end

    #@show hit, miss, new, unbounded
    return verts
end

deleteat(sig, i) = deleteat!(copy(sig), i)


""" generate a random ray orthogonal to the subspace spanned by the given points """
function randray(xs::Points)
    k = length(xs)
    d = length(xs[1])
    v = similar(xs, k-1)

    # Gram Schmidt
    for i in 1:k-1
        v[i] = xs[i] .- xs[k]
        for j in 1:(i-1)
            v[i] = v[i] .- dot(v[i], v[j]) .* v[j]
        end
        v[i] = normalize(v[i])
    end

    u = randn(d)
    for i in 1:k-1
        u = u - dot(u, v[i]) * v[i]
    end
    u = normalize(u)
    return u
end
