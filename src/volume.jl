

""" given two generators `g1`, `g2`, `vertices` of their common boundary
(in generator representation) and the list of all generators, compute the boundary volume.
It works by constructing a polytope from halfspaces between the generators of the boundary
vertices. We use an affine transformation to get rid of the dimension in direction g1->g2
"""
function boundary_area_edges(g1::Int, g2::Int, vertices::AbstractVector{<:Sigma}, generators)
	A = generators[g1]
	B = generators[g2]
	transform = transformation(A, B)
	A = transform(A)
	gen_inds = unique!(sort!(reduce(vcat,vertices)))

	halfspaces = []
	for gen_ind in gen_inds
		gen_ind in (g1,g2) && continue
		B = transform(generators[gen_ind])
		u = normalize(B-A)
		b = dot(u, (A+B)/2)
		hf = HalfSpace(u[2:end], b)
		push!(halfspaces, hf)
	end
	halfspaces = [h for h in halfspaces]

	poly = polyhedron(hrep(halfspaces))
	# performance with libs 6x100
	# nolib: 30 mins
	# CDDLib: 40 mins
	# QHull: 5h (with warnings about not affine polyhedron using a solver)

	vol = hasrays(poly) ? Inf : volume(poly)
	return vol
end

""" affine transformation rotatinig and translating such that the boundary is aligned with
the first dimension. A->B will be mapped to const*[1,0,0,...] and (A+B)/2 to [0,0,...] """
function transformation(A, B)
	R = diagm(ones(length(A)))
	R[:,1] = B-A
	R = inv(qr(R).Q)  # shortcut for gram schmidt orthogonalization
	t = R*(A+B)/2
	transform(x) = R*x - t
end

""" build the connectivity matrix for the SQRA from adjacency and boundary information """
function area_volume(vertices, P::AbstractVector)
	dim = length(P[1])
	conns = adjacency(vertices)

	I = Int[]
	J = Int[]
	V = Float64[]
	Vs = zeros(length(P))

	@showprogress 1 "Voronoi adjacency " for ((g1,g2), sigs) in conns
		A = boundary_area_edges(g1, g2, sigs, P)
		h = norm(P[g1] - P[g2]) / 2
		v = A * h / dim  # volume computation
		push!(I, g1)
		push!(J, g2)
		push!(V, A/h)
		Vs[g1] += v
		Vs[g2] += v
	end
	A = sparse(I, J, V, length(P), length(P))
	A = A + A'

	return A, Vs
end

# for the sqra
function connectivity_matrix(vertices, P::AbstractVector)
	A, Vs = area_volume(vertices, P)
	Vs = replace(Vs, 0 => Inf)  # if we have 0 volume, we have no rates
	Vsi = 1 ./ Vs # TODO: check if we want row or col
	A = A .* Vsi
	return A, Vs
end


""" given vertices in generator-coordinates,
collect the verts belonging to generator pairs, i.e. boundary vertices """
function adjacency(v::Vertices)
	conns = Dict{Tuple{Int,Int}, Vector{Vector{Int}}}()
	#conns=Dict()
	for (sig, r) in v
		for a in sig
			for b in sig
				a >= b && continue
				v = get!(conns, (a,b), [])
				push!(v, sig)
			end
		end
	end
	conns
end

using Polyhedra


""" similar to `boundary_area_edges`, however uses a vector representation and is slower """
function boundary_vrep(g1::Int, g2::Int, inds::AbstractVector{<:Sigma}, vertices::Vertices, generators)
	A = generators[g1]
	B = generators[g2]
	dim = length(A)
	vertex_coords = map(i->vertices[i], inds)
	push!(vertex_coords, A)  # append one voronoi center for full volume
	V = try
			volume(polyhedron(vrep(vertex_coords), QHull.Library()))
		catch e #::QHull.PyCall.PyError
			0
		end

	h = norm(B - A) / 2
	A = dim * V / h
	A, h, V
end
