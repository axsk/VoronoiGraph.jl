@userplot VoronoiPlot

@recipe function f(p::VoronoiPlot; fill_z=nothing)
	v, P = p.args
	fs = faces(v, P)
    label --> ""

	for (i, f) in enumerate(fs)
		@series begin
            if !isnothing(fill_z)
                fill_z := fill_z[i]
            end
			if !hasallrays(f)
				f
			else
				[]
			end
		end
	end
end

function faces(vertices, P::AbstractVector)
	dim = length(P[1])
	conns = Dict{Int, Set{Int}}()
	for (sig, _) in vertices
		for x in sig
			v = get!(conns, x, Set{Int}())
            union!(v, sig)
		end
	end

	faces = [face(i, conns[i], P) for i in 1:length(P)]
	#conns
end

function face(i, js, P)
	halfspaces = []
	A = P[i]
	for j in js
        i == j && continue
		B = P[j]
		u = B-A
		b = dot(u, (A+B)/2)
		hf = HalfSpace(u, b)
		push!(halfspaces, hf)
	end
	halfspaces = [h for h in halfspaces]

	return polyhedron(hrep(halfspaces))
end
