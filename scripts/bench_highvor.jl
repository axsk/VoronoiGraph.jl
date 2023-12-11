using BenchmarkTools
using HighVoronoi
using VoronoiGraph
using MiniQhull
using Plots

function paperplot()
  b = benchsuite(9)
  r = @time run(b)
  try
    plotperf(r)
    savefig("bench.png")
  catch err
    println("plotting failed with $err")
  end
  return r
end

function fastbench()
  run(benchsuite(7), samples=10)
end

function benchsuite(limit=9)
  suite = BenchmarkGroup()

  f(x) = [sum(x .^ 2)]
  f1(x) = sum(x .^ 2)

  for dim in [2, 4, 6]
    for n in [100, 300, 1_000, 3_000, 10_000, 30_000, 100_000, 300_000]
      log10(n) + dim > limit && continue
      local xs = HighVoronoi.VoronoiNodes(rand(dim, n))
      suite[(dim, n, :HighVoronoi)] = @benchmarkable highvor($xs, integrate=false)
      suite[(dim, n, :VoronoiGraph)] = @benchmarkable VoronoiGraph.voronoi($xs)
      suite[(dim, n, :qHull)] = @benchmarkable delaunay($xs)
      #suite[(dim, n,:HighVoronoiMC)] = @benchmarkable VoronoiGeometry($xs, integrator=VI_MONTECARLO, integrand=$f, mc_accurate=(10, 10, 20), silence=true)
      #suite[(dim, n,:VoronoiGraphMC)] = @benchmarkable begin
      #  v = voronoi($xs)
      #  mc_integrate($f1, v[1], $xs, 10, 10)
      #end
    end
  end
  return suite
end

function highvor(xs; kwargs...)
  try
    HighVoronoi.VoronoiGeometry(xs, HighVoronoi.Boundary(); silence=true, kwargs...)
  catch err
    # still getting errors on HighVoronoi@1.2
    println("HV errored with $err")
  end
end

function plotperf(data_dict)
  dims = Set(k[1] for k in keys(data_dict))
  methods = Set(k[3] for k in keys(data_dict))
  times(x) = median(x).time / 1e9
  p = plot()
  linecolor = 1
  for method in methods
    label = method
    for dim in dims
      t = [times(val) for (k, val) in data_dict if k[1] == dim && k[3] == method]
      n = [k[2] for (k, _) in data_dict if k[1] == dim && k[3] == method]
      p = sortperm(n)
      plot!(n[p], t[p]; label, linecolor, marker=:x, markercolor=linecolor)
      label = ""
    end
    linecolor += 1
  end
  plot!(legend=:bottomright, xaxis=:log, yaxis=:log)
end
