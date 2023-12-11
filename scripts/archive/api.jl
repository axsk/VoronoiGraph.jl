

verts(y) = y[1]

vert_sides(y) = reduce(hcat, collect(keys(verts(y)))) :: Matrix

sortverts(x) = reduce(hcat,sort(sort.(eachcol(x))))

function cayleymenger(x::AbstractVector{<:AbstractVector})
    n = length(x)
    M = zeros(n+1, n+1)
    M[1 ,2:end] .= 1
    M[2:end, 1] .= 1
    for i in 1:n, j in (i+1):n
        M[i+1,j+1] = M[j+1,i+1] = norm(x[i]-x[j])^2
    end
    return M
end

# formulas from https://math.stackexchange.com/questions/4056099/circumcenter-of-the-n-simplex/4187480#4187480
function circumcenter(x::AbstractVector{<:AbstractVector})
    M = cayleymenger(x)
    Q = -2 * inv(M)
    q = Q[1,2:end]
    sum(q .* x) / sum(q)
end


function circumcenterdeviations(v,p)
    [norm(circumcenter(p[sig]) - vert) for (sig,vert) in v]
end

function compare_delaunay(x)
    v,p = voronoi(x)
    d = delaunay(x)

    sv = Set(keys(v))
    sq = Set(sort.(eachcol(d)))
    @show length(sv), length(sq)
    @show length(intersect(sv, sq)) / length(sv)
end
