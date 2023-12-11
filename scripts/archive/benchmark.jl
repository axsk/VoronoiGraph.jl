# this was supposed to be run by CI allowing for visual online representation

using BenchmarkTools
using VoronoiGraph
using Random

function bsuite(
    dims=2:5,
    ns=[100, 1000])
    suite = BenchmarkGroup()
    for d in dims
        suite[d] = BenchmarkGroup()
        for n in ns
            Random.seed!(1)
            x = rand(d, n)
            suite[d][n] = @benchmarkable voronoi($x) seconds = 1
        end
    end
    return suite
end

const parameterfile = "../cache/params.json"

function autotune!(suite, load=true)
    if load && isfile(parameterfile)
        try
            println("Loading benchmark parameters")
            loadparams!(suite, BenchmarkTools.load(parameterfile)[1], :evals, :samples)
            return
        catch
        end
    end
    println("Tuning benchmark parameters")
    tune!(suite)
    BenchmarkTools.save("params.json", params(suite))
end

function autorun()
    s = bsuite()
    autotune!(s)
    println("Running benchmark")
    results = run(s, verbose=true)
    println("Saving benchmark results")
    BenchmarkTools.save("output.json", median(results))
    return s, results
end
