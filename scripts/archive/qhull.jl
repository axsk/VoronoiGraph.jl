# performance comparisons to qhull
using BenchmarkTools
using VoronoiGraph
using MiniQhull
using Plots

D = 6
# burn in
for d in 1:D
    voronoi(rand(d, D+1))
end

# times new vs old (pre last 2 commits)
if false
    @benchmark voronoi(rand(2,1000)) # 7.74ms   22.7
    @benchmark voronoi(rand(2,10000)) # 108ms   219
    @benchmark voronoi(rand(5,1000)) # 2.09s    4.91
    @benchmark voronoi(rand(5,2000)) # 5.56s    9.49
    @benchmark voronoi(rand(5,4000)) #         23.6
    @benchmark voronoi(rand(8,100)) # 5.49s    24.5
end



function runmax(f, secs, fallback)
    t = Task(f)
   @async begin
        sleep(secs)
        if !istaskdone(t)
            Base.throwto(t, InterruptException)
        end
    end
    schedule(t)
    try
        fetch(t)
    catch e
        e
        if isa(e, TaskFailedException)
            return fallback
        else
            rethrow(e)
        end
    end
end

function elapsedmax(f, secs)
    t = @elapsed r = runmax(f, secs, -1)
    if r == -1
        return NaN
    else
        return t
    end
end

tmax = 10
plot()
for d in 1:D
    x = rand(d, 1000)
    e1 = elapsedmax(()->delaunay(x), tmax)
    e2 = elapsedmax(()->voronoi(x), tmax)
    scatter!([d], [e1 e2], labels=["D" "V"])
end

plot!()
