const Point{T} = AbstractVector{T} where T<:Real
const Points = AbstractVector{<:Point}
const Sigma = AbstractVector{<:Integer}  # Encoded the vertex by the ids of its generators.
const Vertex = Tuple{<:Sigma, <: Point}
const Vertices = Dict{<:Sigma, <:Point}

dim(xs::Points) = length(xs[1])

voronoi(x; kwargs...) = voronoi(vecvec(x); kwargs...)

""" construct the voronoi diagram from `x` through breadth-first search """
function voronoi(xs::Points, searcher = Raycast(xs))
    sig, r = descent(xs, searcher)
    verts, rays = explore(sig, r, xs, searcher)
    return verts::Vertices, xs, rays
end

voronoi_random(x, args...; kwargs...) = voronoi_random(vecvec(x), args...; kwargs...)

""" construct a (partial) voronoi diagram from `x` through a random walk """
function voronoi_random(xs::Points, iter=1000; maxstuck=typemax(Int))
    searcher = Raycast(xs)
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

    sig = SVector{dim(xs)+1}(sig)
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
        randdir = rand(1:length(sig))
        sig, r = walkray(sig, r, xs, searcher, randdir)
        r == Inf && continue  # we encountered an infinite ray
        get!(verts, sig) do
            nonew = 0
            return r
        end
        nonew > maxstuck && break
    end

    return verts
end


""" find the vertex connected to `v` by moving away from its `i`-th generator """
function walkray(sig::Sigma, r::Point, xs::Points, searcher, i)
    sig_del = deleteat(sig, i)
    u = u_default(sig, xs, i)
    sig′, t = raycast(sig_del, r, u, xs, searcher)
    if t < Inf
        r′ = r + t*u
        return sig′, r′
    else
        return sig, r  # if the vertex has an unbounded ray, return the same vertex
    end
end



function u_compare(sig, xs, i)
    u1 = u_randray(sig, xs, i)
    u2 = u_qr(sig, xs, i)
    #u3 = u_qr_short(deleteat(sig, i), xs, sig[i])
    u3=u2
    if !(isapprox(u1, u2) && isapprox(u1, u2))
        @warn "Not approx" u1 u2 u3
    end
    return u1
end

function u_randray(sig, xs, i)
    sig_del = deleteat(sig, i)
    u = randray_modified(xs[sig_del])
    if (u' * (xs[sig[i]] - xs[sig_del[1]])) > 0
        u = -u
    end
    return u
end

# replace gram-schmidt of randray by qr
function u_qr(sig, xs, i)
    n = length(sig)
    d = length(xs[1])
    X = MMatrix{d, n-1, eltype(xs[1])}(undef)
    for j in 1:i-1
        X[:, j] = xs[sig[j]]
    end
    for j in i:n-1
        X[:, j] = xs[sig[j+1]]
    end
    origin = X[:, end]
    X[:, end] = xs[sig[i]]
    X .-= origin
    X = SMatrix(X)
    Q, R = qr(X)
    u = -Q[:,end]  * sign(R[end,end])
    return u
end

# same as u_qr, shorter and slower
function u_qr_short(sig_del, xs, opposing)
    sig = vcat(sig_del[2:end], [opposing])
    origin = xs[sig_del[1]]
    X = mapreduce(i->xs[i] - origin,hcat, sig)
    X = Array(X)
    #@show t
    q = qr(X)
    u = q.Q[:,end]
    u *= -sign(q.R[end,end])

    return u
end

global hit::Int
global miss::Int
global new::Int
global cache::Bool = true

""" BFS of vertices starting from `S0` """
function explore(sig, r, xs::Points, searcher) # :: Vertices
    verts = Dict(sig=>r)
    queue = [sig=>r]
    edgecount = Dict{SVector{dim(xs), Int64}, Int}()
    sizehint!(edgecount, @show round(Int, expected_edges(xs)))
    sizehint!(verts, round(Int, expected_vertices(xs)))
    rays = Pair{Vector{Int64}, Int}[]

    global hit = 0
    global miss = 0
    global new = 0
    unbounded = 0

    while length(queue) > 0
        #yield() we used this to ease ctrl+c, but it takes ~10% runtime
        (sig,r) = pop!(queue)
        for i in 1:length(sig)

            if cache && get(edgecount, deleteat(sig, i), 0) == 2
                hit += 1
                continue
            end

            sig′, r′ = walkray(sig, r, xs, searcher, i)

            if sig′ == sig
                push!(rays, sig => i)
                unbounded += 1
                continue
            end

            if !haskey(verts, sig′)
                push!(queue, sig′ => r′)
                push!(verts, sig′ => r′)
                if cache
                    for j in 1:length(sig′)
                        edge = deleteat(sig′, j)
                        edgecount[edge] = get(edgecount, edge, 0) + 1
                    end
                end
                new += 1
            else
                miss += 1
            end
        end
    end
    @show length(edgecount)
    @show Base.summarysize(verts) / 1024^2
    @show Base.summarysize(edgecount) / 1024^2

    @show hit, miss, new, unbounded
    return verts, rays
end

const u_default = u_qr

deleteat(sig::Vector, i) = deleteat!(copy(sig), i)
deleteat(x,y) = StaticArrays.deleteat(x,y)

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

""" generate a random ray orthogonal to the subspace spanned by the given points """
function randray_modified(xs::Points)
    k = length(xs)
    v = collect(xs)
    for i in 1:k
        v[i] -= v[k] # subtract origin
    end
    v[k] = randn(dim(xs))

    # modified Gram-Schmidt (for better stability)
    for i in 1:k
        v[i] = normalize(v[i])
        for j in i+1:k
            v[j] -= dot(v[i], v[j]) * v[i]
        end
    end

    u = v[k]
    return u
end

"""
vertexheuristic(d, n)

expected number of vertices for `n` points in `d` dimensions
c.f. [Dwyer, The expected number of k-faces of a Voronoi Diagram (1993)] """
function expected_vertices(d, n)

    # from the previous work, only good for d>5
    # 2^((d+1)/2) * exp(1/4) * pi^((d-1)/2) * d^((d-2)/2) / (d+1) * n

    # lower and upperbound coincide for k==d
    lowerbound(d, d) * n
end

function expected_edges(d, n)
    expected_vertices(d, n) * (d+1) / 2
end

expected_edges(xs::Points) = expected_edges(dim(xs), length(xs))
expected_vertices(xs::Points) = expected_vertices(dim(xs), length(xs))

# Theorem 1
# For fixed k and d as n grow without bound, EF(k,d)  ̃ C(k,d)n with C(k,d) ≥ lowerbound(k,d)
function lowerbound(k, d)
    2 * pi^(k/2) * d^(k-1) / (k*(k+1)) * beta(d*k/2, (d-k+1)/2) ^ (-1) * (gamma(d/2) / gamma((d+1)/2))^k
end

# Theorem 1
# For fixed k and d as n grow without bound, EF(k,d)  ̃ C(k,d)n with C(k,d) ≤ upperbound(k,d)
function upperbound(k, d)
    2 * d^(d-2) * pi^((d-1)/2) * gamma((d^2+1)/2) * gamma(d/2)^d /
    ((d+1) * (d+2) * gamma(d^2/2) * gamma((d+1)/2)^d) * binomial(d+2, k+1)
end
