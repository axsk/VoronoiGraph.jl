
## Implementations of different raycasting search algorithms

# default raycast
Raycast(xs) = RaycastCompare(KDTree(xs), 100000, 1e-8, zeros(4))

struct SearchBruteforce
    tmin
end

""" shooting a ray in the given direction, find the next connecting point.
This is the bruteforce variant, using a linear search to find the closest point """
function raycast(sig::Sigma, r, u, xs, searcher::SearchBruteforce)
    (tau, ts) = [0; sig], Inf
    x0 = xs[sig[1]]

    for i in 1:length(xs)
        i in sig && continue
        x = xs[i]
        t = (sum(abs2, r .- x) - sum(abs2, r .- x0)) / (2 * u' * (x-x0))
        if searcher.tmin < t < ts
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
        return [0; sig], Inf
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


struct RaycastCompare
    tree::KDTree
    tmax
    eps
    timings
end

function raycast(sig::Sigma, r::Point, u::Point, xs::Points, searcher::RaycastCompare)
    s1 = SearchBruteforce(searcher.eps)
    s2 = SearchBisection(searcher.tree, searcher.tmax, searcher.eps)
    s3 = SearchIncircle(searcher.tree, searcher.tmax)
    s4 = SearchIncircleSkip(searcher.tree)

    t1 = @elapsed r1 = raycast(sig, r, u, xs, s1)
    t2 = @elapsed r2 = raycast(sig, r, u, xs, s2)
    t3 = @elapsed r3 = raycast(sig, r, u, xs, s3)
    t4 = @elapsed r4 = raycast(sig, r, u, xs, s4)

    searcher.timings .+= [t1, t2, t3, t4]

    if !(r1[1]==r2[1]==r3[1]==r4[1])
        @warn "raycast algorithms return different results" r1 r2 r3 r4
    end
    return r4
end