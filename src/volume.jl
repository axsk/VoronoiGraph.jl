function boundary_area(args...)
    error("This function requires the Polyhedron.jl package. Please install it and try again.")
end

""" affine transformation rotatinig and translating such that the boundary is aligned with
the first dimension. A->B will be mapped to const*[1,0,0,...] and (A+B)/2 to [0,0,...] """
function transformation(A, B)
    R = diagm(ones(length(A)))
    R[:,1] = B-A
    R = collect(inv(qr(R).Q))  # Shortcut for Gram-Schmidt orthogonalization.
    t = R*(A+B)/2
    transform(x) = R*x - t
    return transform
end

function transformation(A::SVector{N, T}, B::SVector{N, T}) where {N, T}
    R = MMatrix{N,N,T}(I)
    R[:,1] = B-A
    R = inv(qr(R).Q)
    t = R * (A+B)/2
    transform(x) = R*x - t
end

""" build the connectivity matrix for the SQRA from adjacency and boundary information """
function volumes(vertices, P::AbstractVector)
    dim = length(P[1])
    conns = adjacency(vertices)

    I = Int[]
    J = Int[]
    As = Float64[]
    Vs = zeros(length(P))

    #=@showprogress 1 "Voronoi volumes "=# for ((g1,g2), sigs) in conns
        a = boundary_area(g1, g2, sigs, P)
        h = norm(P[g1] - P[g2]) / 2
        v = a * h / dim  # Volume computation
        push!(I, g1)
        push!(J, g2)
        push!(As, a)
        Vs[g1] += v
        Vs[g2] += v
    end
    A = sparse(I, J, As, length(P), length(P))
    A = A + A'

    return A, Vs
end

""" given vertices in generator-coordinates,
collect the verts belonging to generator pairs, i.e. boundary vertices """
function adjacency(v::Vertices)
    conns = Dict{Tuple{Int,Int}, Vector{Vector{Int}}}()

    for (sig, r) in v
        for a in sig
            for b in sig
                a >= b && continue
                v = get!(conns, (a,b), [])
                push!(v, sig)
            end
        end
    end

    return conns
end
