
## Implementations of different raycasting search algorithms

# default raycast
Raycast(xs) = RaycastIncircleSkip(KDTree(xs))

struct RaycastBruteforce end

""" shooting a ray in the given direction, find the next connecting point.
This is the bruteforce variant, using a linear search to find the closest point """
function raycast(sig::Sigma, r, u, xs, searcher::RaycastBruteforce)
    (tau, ts) = [0; sig], Inf
    x0 = xs[sig[1]]

    c = maximum(dot(xs[g], u) for g in sig)
    skip(i) = (dot(xs[i], u) <= c) || i ∈ sig


    for i in 1:length(xs)
        skip(i) && continue
        x = xs[i]
        t = (sum(abs2, r .- x) - sum(abs2, r .- x0)) / (2 * u' * (x-x0))
        if 0 < t < ts
            (tau, ts) = vcat(sig, [i]), t
        end
    end

    return sort(tau), ts
end


struct RaycastBisection{T}
    tree::T
    tmax
    eps::Float64
end

""" shooting a ray in the given direction, find the next connecting point.
This variant (by Poliaski, Pokorny) uses a binary search """
function raycast(sig::Sigma, r::Point, u::Point, xs::Points, searcher::RaycastBisection)
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


struct RaycastIncircle{T}
    tree::T
    tmax::Float64
end

""" Shooting a ray in the given direction, find the next connecting point.
This variant uses an iterative NN search """
function raycast(sig::Sigma, r::Point, u::Point, xs::Points, searcher::RaycastIncircle)
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


struct RaycastIncircleSkip{T}
    tree::T
end

"""
raycast(sig::Sigma, r::Point, u::Point, xs::Points, seacher::RaycastIncircleSkip)

Return `(tau, t)` where `tau = [sig..., i]` are the new generators
and `t` the distance to walk from `r` to find the new representative.

`::RaycastIncircleSkip` uses the nearest-neighbours skip predicate
to find the initial candidate on the right half plane.
Resumes with successive incircle searches until convergence
(just as in ``::RaycastIncircle`).

Q: Why is this routine split in these two parts
instead of only iterating on the right half plane?
A: The first might find no generator in case of a ray. This is handled explicitly here.
"""

function raycast(sig::Sigma, r::Point, u::Point, xs::Points, searcher::RaycastIncircleSkip)

    x0 = xs[sig[1]]

    # only consider points on the right side of the hyperplane
    c = maximum(dot(xs[g], u) for g in sig)
    skip(i) = (dot(xs[i], u) <= c) || i ∈ sig

    # shift candidate onto the plane spanned by the generators
    candidate = r + u * (u' * (x0-r))

    # compute heuristic t assuming the resulting delauney simplex was regular
    # this reduces number of extra searches by about 10%
    if length(sig) > 1
        n = length(sig)
        radius = norm(candidate-x0)
        #radius = sum(norm(candidate-xs[s]) for s in sig) / n  # better but slower
        t = radius / sqrt((n+1)*(n-1))
        candidate += t * u
    end

    is, _ = knn(searcher.tree, candidate, 1, false, skip)
    (length(is) == 0) && return [0; sig], Inf  # no point was found
    i = is[1]

    # sucessively reduce incircles unless nothing new is found
    while true
        x = xs[i]
        t = (sum(abs2, r - x) - sum(abs2, r - x0)) / (2 * u' * (x-x0))
        candidate = r + t*u
        j, d = nn(searcher.tree, candidate)

        (j in sig || j == i) && break  # converged to the smallest circle

        dold = sqrt(sum(abs2, x0-candidate))
        isapprox(d, dold) && @warn "degenerate vertex at $sig + [$i] ($d $dold)"

        i = j
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

RaycastCompare(xs) = RaycastCompare(KDTree(xs), 1_000., 1e-8, zeros(4))

function raycast(sig::Sigma, r::Point, u::Point, xs::Points, searcher::RaycastCompare)
    s1 = RaycastBruteforce()
    s2 = RaycastBisection(searcher.tree, searcher.tmax, searcher.eps)
    s3 = RaycastIncircle(searcher.tree, searcher.tmax)
    s4 = RaycastIncircleSkip(searcher.tree)

    t1 = @elapsed r1 = raycast(sig, r, u, xs, s1)
    t2 = @elapsed r2 = raycast(sig, r, u, xs, s2)
    t3 = @elapsed r3 = raycast(sig, r, u, xs, s3)
    t4 = @elapsed r4 = raycast(sig, r, u, xs, s4)

    searcher.timings .+= [t1, t2, t3, t4]

    if !(r1[1]==r2[1]==r3[1]==r4[1]) && r3[2] < searcher.tmax
        @warn "raycast algorithms return different results" r1 r2 r3 r4 tuple(r...)
    end
    return r4
end
