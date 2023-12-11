using VoronoiGraph
using Plots

d=1:16
c=[VoronoiGraph.lowerbound(d,d) for d in d]

plot(d, c, label="", yaxis=:log)
