### Count the number of nearest neighbour calls

using VoronoiGraph
using NearestNeighbors

time_from_now(seconds) = round(Int, 10^9 * seconds + time_ns())

## Container for a nearest neighbour search counting the number of nn calls
mutable struct TimedCountingTree{T}
    tree::T
    count::Int
    killtime::Int
end

function NearestNeighbors.knn(ct::TimedCountingTree, args...)
    if time_ns() > ct.killtime
        error("time limit exceeded")
    end
    ct.count += 1
    return NearestNeighbors.knn(ct.tree, args...)
end

function NearestNeighbors.nn(ct::TimedCountingTree, point)
    if time_ns() > ct.killtime
        error("time limit exceeded")
    end
    ct.count += 1
    return NearestNeighbors.nn(ct.tree, point)
end


# comparison of the two calls

function countnn_ris(xs; maxtime=60, kwargs...)
    s = TimedCountingTree(KDTree(xs), 0, time_from_now(maxtime))
    local v
    try
        v, = voronoi(VoronoiGraph.vecvec(xs), VoronoiGraph.RaycastIncircleSkip(s))
    catch
        return NaN
    end
    return s.count / @show length(v)
end

function countnn_bis(xs; tmax=10, eps=1e-3, maxtime=60, kwargs...)
    s = TimedCountingTree(KDTree(xs), 0, time_from_now(maxtime))
    local v
    try
        v, = voronoi(VoronoiGraph.vecvec(xs), VoronoiGraph.RaycastBisection(s, tmax, eps))
    catch
        return NaN
    end
    return s.count / @show length(v)
end

using Random

function tabular_comparison(dims=[2, 3, 4, 5, 6, 7], n=1000, maxtime=300, seed=3)
    Random.seed!(seed)
    calls = zeros(length(dims), 4)
    println("maximal time: $(length(calls)*maxtime) seconds")
    for (i, dim) in enumerate(dims)
        xs = rand(dim, n)
        @time calls[i, 1], = countnn_bis(xs, eps=1e-4; maxtime)
        @time calls[i, 2], = countnn_bis(xs, eps=1e-8; maxtime)
        VoronoiGraph.use_heuristic(false)
        @time calls[i, 3], = countnn_ris(xs; maxtime)
        VoronoiGraph.use_heuristic(true)
        @time calls[i, 4], = countnn_ris(xs; maxtime)
    end
    return calls
end
using Printf
function print_latex(x, captions=["bisection (\$1e-4\$)", "bisection (\$1e-8\$)", "incircle", "incircle heuristic"])
    d = range(2, size(x, 2) + 1)
    println("\\begin{tabular}{l|c|c|c|c|}")
    print("method & " * join(["\$d=$i\$" for i in d], " & ") * "\\\\ \\hline\n")

    for (c, r) in zip(captions, eachrow(x))
        print(c)
        for x in r
            @printf("& %.2f ", x)
        end
        println("\\\\")
    end
    println("\\end{tabular}")
end
