using HighVoronoi
using Plots

function integrand(periods)
    fevals = [0]
    function f(x)
        fevals .+= 1
        [sin(sum(x .^ 2) / dim * periods * 2 * 2pi)]
    end
    return f, fevals
end

function voronoi(xs; integrator, nmc=100, periods=1)
    f, fevals = integrand(periods)
    v = @time VoronoiGeometry(VoronoiNodes(xs), cuboid(size(xs, 1)), integrator=integrator, integrand=f, mc_accurate=(nmc, 3, 2), silence=true)
    nverts = length(reduce(union, keys.(v.Integrator.Integral.MESH.All_Verteces)))
    @show nverts, nverts / size(xs, 2)
    @show fevals
    v
end

montecarlo(xs, nmc=100) = voronoi(xs, integrator=VI_MONTECARLO, nmc=nmc)
heuristic(xs) = voronoi(xs, integrator=VI_HEURISTIC)
polygon(xs) = voronoi(xs, integrator=VI_POLYGON)


function compare_mc_h(;
    dim=3,
    n=1000,
    xs=rand(dim, n),
    nmc=100,
    nmc2=1000,
    periods=10,
    obs=intsurf)

    #xs = rand(dim, n)
    mc = montecarlo(xs, nmc)
    mc2 = montecarlo(xs, nmc2)
    h = heuristic(xs)

    plot()
    plot_rel_errors(mc2, mc, obs=obs, label="mc")
    plot_rel_errors(mc2, h, obs=obs, label="h")

    (; xs, mc, mc2, h)

end

function l2loss(truth, est; obs)
    y = obs(truth)
    ŷ = obs(est)
    inds = map(!isnan, y + ŷ)
    sqrt(mean((y[inds] .- ŷ[inds]) .^ 2))
end

using StatsBase
function plot_rel_errors(truth, est; obs=area, nbins=20, q=0.1, abs=false, label="", scatter=false)
    truth = obs(truth)
    if abs
        err = obs(est) - truth
    else
        err = obs(est) ./ truth
    end
    lims = quantile(filter(isfinite, truth), [q, 1 - q])
    bins = range(lims..., length=nbins)

    means = zeros(length(bins) - 1)
    stds = zeros(length(bins) - 1)
    for i in 1:length(bins)-1
        j = findall(truth .> bins[i] .&& truth .< bins[i+1])
        e = filter(isfinite, err[j])
        means[i] = mean(filter(isfinite, err[j]))
        stds[i] = sqrt(mean((e .- 1) .^ 2))
    end
    scatter && scatter!(truth, rel, xlims=lims)
    plot!(bins[1:end-1] .+ step(bins), means, label="$label mean")
    plot!(bins[1:end-1] .+ step(bins), stds, label="$label std")
    ylims!(0, 1.2) |> display

    return means, stds
end

function plot3h()
    xs = rand(3, 10000)
    @time r = [compare_mc_h(xs=xs, nmc=10_000, nmc2=100_000, periods=2) for i in 1:2]
end



area(v::VoronoiGeometry) = v.Integrator.Integral.area |> sumundef
vol(v::VoronoiGeometry) = v.Integrator.Integral.volumes |> sumundef
intsurf(v::VoronoiGeometry) = v.Integrator.Integral.interface_integral |> sumundefsurf
intvol(v::VoronoiGeometry) = v.Integrator.Integral.bulk_integral |> sumundef

function sumundef(a)
    a = copy(a)
    for i in eachindex(a)
        isassigned(a, i) && continue
        a[i] = [NaN]
    end
    sum.(a)
end

function sumundefsurf(a)
    a = copy(a)
    for i in eachindex(a)
        isassigned(a, i) && continue
        a[i] = [[NaN]]
    end
    map(x -> sum(x[1]), a)
end
