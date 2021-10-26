const Point{T} = AbstractVector{T} where T<:Real
const Points = AbstractVector{<:Point}
const Sigma = AbstractVector{<:Integer}  # Encoded the vertex by the ids of its generators.
const Vertex = Tuple{<:Sigma, <: Point}
const Vertices = Dict{<:Sigma, <:Point}


voronoi(x; kwargs...) = voronoi(vecvec(x); kwargs...)

""" construct the voronoi diagram from `x` through breadth-first search """
function voronoi(xs::Points)
    #searcher = SearchIncircle(tmax, KDTree(xs))
    searcher = RaycastCompare(KDTree(xs), 10000, 1e-8, zeros(4))
    sig, r = descent(xs, searcher)
    verts = explore(sig, r, xs, searcher)
    @show searcher.timings
    return verts::Vertices, xs
end

voronoi_random(x, args...; kwargs...) = voronoi_random(vecvec(x), args...; kwargs...)

""" construct a (partial) voronoi diagram from `x` through a random walk """
function voronoi_random(xs::Points, iter=1000; tmax=1000, maxstuck=typemax(Int))
    searcher = SearchIncircle(KDTree(xs), tmax)
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


## Implementations of different raycasting search algorithms

Raycast(xs) = RaycastCompare

struct SearchBruteforce
    eps
end

""" shooting a ray in the given direction, find the next connecting point.
This is the bruteforce variant, using a linear search to find the closest point """
function raycast(sig::Sigma, r, u, xs, searcher::SearchBruteforce)
    eps = searcher.eps
    (tau, ts) = [0; sig], Inf
    x0 = xs[sig[1]]

    for i in 1:length(xs)
        i in sig && continue
        x = xs[i]
        t = (sum(abs2, r .- x) - sum(abs2, r .- x0)) / (2 * u' * (x-x0))
        if eps < t < ts
            (tau, ts) = vcat(sig, [i]), t
        end
    end

    return sort(tau), ts
end


struct SearchBisection
    tree::KDTree
    tmax
    eps::Float64
end

""" shooting a ray in the given direction, find the next connecting point.
This variant (by Poliaski, Pokorny) uses a binary search """
function raycast(sig::Sigma, r::Point, u::Point, xs::Points, searcher::SearchBisection)
    tau, tl, tr = [], 0, searcher.tmax
    x0 = xs[sig[1]]

    while tr-tl > searcher.eps
        tm = (tl+tr)/2
        i, _ = nn(searcher.tree, r+tm*u)
        x = xs[i]
        if i in sig
            tl = tm
        else
            tr = (sum(abs2, r .- x) - sum(abs2, r .- x0)) / (2 * u' * (x-x0))
            tau = vcat(sig, [i])

            # early stopping
            idxs, _ = knn(searcher.tree, r+tr*u, length(sig)+1, true)
            length(intersect(idxs, [sig; i])) == length(sig)+1 && break
        end
    end

    if tau == []
        tau = [0; sig]
        tr = Inf
    end

    return sort(tau), tr
end


struct SearchIncircle
    tree::KDTree
    tmax::Float64
end

""" Shooting a ray in the given direction, find the next connecting point.
This variant uses an iterative NN search """
function raycast(sig::Sigma, r::Point, u::Point, xs::Points, searcher::SearchIncircle)
    i = 0
    t = 1
    x0 = xs[sig[1]]
    local d, n

    # find a t large enough to include a non-boundary (sig) point
    while t < searcher.tmax
        n, d = nn(searcher.tree, r+t*u)
        if d==Inf
            warn("d==Inf in raycast expansion, this should never happen")
            return [0; sig], Inf
        end

        if n in sig
            t = t * 2
        else
            i = n
            break
        end
    end

    if i == 0
        return [0; sig], Inf
    end

    # sucessively reduce incircles unless nothing new is found
    while true
        x = xs[i]
        t = (sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        j, _ = nn(searcher.tree, r+t*u)
        if j in [sig; i]
            break
        else
            i = j
        end
    end

    tau = sort([i; sig])

    return tau, t
end

struct RaycastCompare
    tree::KDTree
    tmax
    eps
    timings
end

function raycast(sig::Sigma, r::Point, u::Point, xs::Points, searcher::RaycastCompare)
    s1 = SearchIncircle(searcher.tree, searcher.tmax)
    s2 = SearchIncircleSkip(searcher.tree)
    s3 = SearchBisection(searcher.tree, searcher.tmax, searcher.eps)
    s4 = SearchBruteforce(searcher.eps)

    t1 = @elapsed r1 = raycast(sig, r, u, xs, s1)
    t2 = @elapsed r2 = raycast(sig, r, u, xs, s2)
    t3 = @elapsed r3 = raycast(sig, r, u, xs, s3)
    t4 = @elapsed r4 = raycast(sig, r, u, xs, s4)

    searcher.timings .+= [t1, t2, t3, t4]

    #@show r1, r2
    if !(r1[1]==r2[1]==r3[1]==r4[1])
        @show r1, r2, r3, r4
        error("raycast algoriths return different results")
    end
    #@assert r1[1] == r2[1]
    return r1
end


#=
function raycast_compare(sig, r, u, P, searcher)
    r1  = raycast_incircle(sig,r,u,P,searcher)
    r2  = raycast_intersect(sig,r,u,P,searcher)
    r3  = raycast_bruteforce(sig,r,u,P)
    @assert r1[1] == r2[1]
    @assert r2[1] == r3[1]
    return r1
end
=#

struct SearchIncircleSkip
    tree::KDTree
end

function raycast(sig::Sigma, r::Point, u::Point, xs::Points, searcher::SearchIncircleSkip)

    x0 = xs[sig[1]]

    #c = dot(x0, u)
    c = maximum(dot(xs[g], u) for g in sig)
    #c = maximum(dists)
    #@show c

    # only consider points on the right side of the hyperplane
    #skip(i) = (u' * (xs[i]-x0) <= 0) || i ∈ sig
    skip(i) = (dot(xs[i], u) <= c) || i ∈ sig

    # try catch workaround for https://github.com/KristofferC/NearestNeighbors.jl/issues/127
    local i, t
    try
        i, t = nn(searcher.tree, r + u * (u' * (x0-r)), skip)
    catch
        return [0, sig], Inf
    end

    t == Inf && return [0; sig], Inf

    # sucessively reduce incircles unless nothing new is found
    while true
        x = xs[i]
        t = (sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        j, _ = nn(searcher.tree, r+t*u)
        if j in [sig; i]
            break
        else
            i = j
        end
    end

    tau = sort([i; sig])

    return tau, t
end
