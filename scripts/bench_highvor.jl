using BenchmarkTools
using HighVoronoi
using VoronoiGraph
using MiniQhull
using Plots

function benchsuite(limit=7)
  suite = BenchmarkGroup()

  f(x) = [sum(x .^ 2)]
  f1(x) = sum(x .^ 2)

  for dim in [2, 6]
    for n in [100, 300, 1_000, 3_000, 10_000, 30_000, 100_000, 300_000]
      log10(n) + dim > limit && continue
      local xs = HighVoronoi.VoronoiNodes(rand(dim, n))
      suite[(dim, n)] = BenchmarkGroup(["($dim, $n)"])

      suite[(dim, n)]["HighVoronoi"] = @benchmarkable HighVoronoi.VoronoiGeometry($xs, Boundary(), silence=true)
      suite[(dim, n)]["VoronoiGraph"] = @benchmarkable VoronoiGraph.voronoi($xs)
      suite[(dim, n)]["qHull"] = @benchmarkable delaunay($xs)

      suite[(dim, n)]["HighVoronoi MC"] = @benchmarkable VoronoiGeometry($xs, Boundary(), integrator=VI_MONTECARLO, integrand=$f, mc_accurate=(10, 10, 20), silence=true)
      suite[(dim, n)]["VoronoiGraph MC"] = @benchmarkable begin
        v = voronoi($xs)
        mc_integrate($f1, v[1], $xs, 10, 10)
      end
    end
  end
  return suite
end

#=

function time_mc(dim=4, n=20000)
  ns = VoronoiNodes(rand(dim, n))
  @time VoronoiGeometry(ns, Boundary(), integrator=VI_MONTECARLO, integrand=f, mc_accurate=(10, 10, 20), silence=true)
  @time begin
    v = voronoi(ns)
    mc_integrate(f1, v[1], ns, 10, 10)
  end
  nothing
end




ratios = map(collect(r)) do ((dim, nx), r)
  (dim, nx) => /(median.([r["HighVoronoi"].times, r["VoronoiGraph"].times])...)
end

sort(collect(r), by=first)

=#
#=
julia> sort(collect(r), by=first)
20-element Vector{Pair{Tuple{Int64, Int64}, Float64}}:
    (2, 100) => 7.016929340643119
   (2, 1000) => 5.421190754929933
  (2, 10000) => 6.007883866403234
  (2, 30000) => 6.6676912093263105
 (2, 100000) => 7.586856294087824
    (3, 100) => 4.28197215146044
   (3, 1000) => 3.682231540026104
  (3, 10000) => 4.011504655703367
  (3, 30000) => 3.96267038467643
 (3, 100000) => 4.545421601253029
    (4, 100) => 3.225445017795933
   (4, 1000) => 3.400324085185257
  (4, 10000) => 3.587482731340586
  (4, 30000) => 3.920902844738698
 (4, 100000) => 4.092963695283695
    (5, 100) => 3.356166405797583
   (5, 1000) => 3.318095967735031
  (5, 10000) => 3.575366420291789
    (6, 100) => 3.247622761266849
   (6, 1000) => 3.3487722004773985
   =#

#methods_int = ["HighVoronoi MC", "VoronoiGraph MC"]

dimkeys(r, d) = sort(collect(keys(filter(r) do pair
  pair[1][1] == d
end)))

using Plots
function plotperf(r; d=2, methods=["qHull", "HighVoronoi", "VoronoiGraph"], kwargs...)
  plot()
  d = 2
  ks = dimkeys(r, d)
  suff = ""

  for (i, method) in enumerate(methods)
    n = []
    ts = []
    i > 1 && (suff = ".jl")
    for k in ks
      push!(n, k[2])
      push!(ts, median(r[k][method]).time / 1e9)
    end
    plot!(n, ts, label="$method$suff"; linecolor=i, kwargs...)
  end

  for d in [6]

    @show ks = dimkeys(r, d)

    for (i, method) in enumerate(methods)
      n = []
      ts = []
      for k in ks
        push!(n, k[2])
        push!(ts, median(r[k][method]).time / 1e9)
      end
      plot!(n, ts, label=""; linecolor=i, kwargs...)
    end
  end


  plot!(xaxis=:log, yaxis=:log, title="Compute time", xlabel="n", ylabel="time (s)")
  annotate!(4.6, -1.5, "d = 2")
  annotate!(2.5, 1.2, "d = 6")
  #savefig("perf.pdf")
end

function runplots(r=BenchmarkTools.run(benchsuite()))
  plotperf(r, marker=:o, markersize=2) |> display
  return r
end

# TO RUN THE BENCHMARKS:
#run(benchsuite())