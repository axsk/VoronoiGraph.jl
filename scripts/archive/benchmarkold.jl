using VoronoiGraph

function benchmark(d, n, mc=0)
    data = rand(d, n)
    t1 = @elapsed v, p = voronoi(data)
    if mc > 0
        t2 = @elapsed V, A = mc_volumes(v, p, mc)
    else
        t2 = @elapsed V, A = volumes(v, p)
    end
    return t1, t2
end

using Plots



function plotbench(; ds=2:5, n=1000, mc=0, nd=0)
    ts = mapreduce(vcat, ds) do d
        n = nd > 0 ? nd^d : n
        t1, t2 = @time benchmark(d, n, mc)
        [t1 t2]
    end
    plot(ds, ts, yaxis = :log) |> display
    return ts
end

lucaplot() = plotbench(ds=2:5, nd=5, mc=0)

quickplot() = plotbench(ds=2:5, n=100, mc=100)
using LinearAlgebra: norm




""" return all vertices having rays """
boundaryverts(rays) = unique(v[1] for v in rays)

""" return all points lying at the boundary,
i.e. whose cell is delimited by a vertex having a (infinite) ray
Note: The corresponding cells can still be bounded, but do neighbour to an unbounded cell
"""
boundarypoints(rays) = reduce((x,y)->unique(vcat(x,y)), v[1] for v in rays)

function innerpoints(p, rays)
    b = boundarypoints(rays)
    inds = setdiff(1:length(p), b)
    return inds
end

function innerverts(v, rays)
    b = boundarypoints(rays)
    vs = filter(s->!all(i -> i in b, s), keys(v))
    return vs
end

function innerrelation(v, p, rays)
    ip = length(innerpoints(p, rays))
    iv = length(innerverts(v, rays))
    return ip, iv, iv/ip
end

function plot_innerrelation(d=1:3, ns=10:10:1000)
    p = plot(title = "inner relation", xlabel="n", ylabel="#v/#p", yaxis=:log)
    for d in d
        ir = []
        br = []
        for n in ns
            @time "d=$d n=$n" begin
                vor = voronoi(rand(d,n) .*2 .- 1)
                push!(ir, innerrelation(vor...)[3])
                push!(br, boxrelation(vor[1], vor[2])[3])
            end
        end

        plot!(p, ns, ir, label="d=$d (inner)", color=d)
        plot!(p, ns, br, label="d=$d (box)", color=d, linestyle=:dash)
    end
    p
end


""" estimate the ratio of vertices inside a smaller box volume """
function boxrelation(v, p, r=.5)
    ip = count(x->norm(x, Inf)<r, p)
    iv = count(x->norm(x, Inf)<r, values(v))
    return ip, iv, iv/ip
end

function vratios(d, n, r=.5)
    v, p = voronoi(rand(d, n))
    scatter!([length(p)], [vertratio(v, p, r)], color=length(p[1]), xaxis=:log)
end

function boxsize(p, rays)
    binds = boundarypoints(rays)
    r = minimum(x->norm(x, Inf), @view p[binds])
end



function mc_variance(f, i, xs, n)
    coll = Array{Any}(undef, n, 4)
    for j in 1:n
        coll[j,:] .= mc_integrate(f, i, xs, 1, 1)
    end
#return coll

    my = sum(coll[:,3]) / n
    vy = sum((coll[:,3] .- my).^2) / n

    mdy, vdy = variance_vec_of_spvec(coll[:,4])

    return my, sqrt.(vy) ./ my, mdy, vdy
end

function variance_vec_of_spvec(vec)
    n = length(vec)
    mn = sum(vec) / n
    var = zero(mn)
    cnt = zero(mn)
    for s in vec
        i = s.nzind[1]
        val = s.nzval[1]
        cnt[i] += 1
        var[i] += (val - mn[i]) ^ 2
    end
    var.nzval ./= cnt.nzval
    return mn, var
end
