""" given a a voronoi diagram (vs,xs,rays)
return the indices of all cells with bounded volume """
function boundedcells(verts, xs, rays)
    unbounded = unique(reduce(vcat, first.(rays)))
    [i for i in 1:length(xs) if !in(i, unbounded)]
end



""" given a voronoi diagram (vs, xs, rays)
return the indices of all cells completely contained inside the unit box """
function containedcells(verts, xs, rays)
    outside = [k for (k, v) in verts if !all(0 .< v .< 1)]
    outside = unique(reduce(vcat, outside))
    [i for i in 1:length(xs) if !in(i, outside)]
end

