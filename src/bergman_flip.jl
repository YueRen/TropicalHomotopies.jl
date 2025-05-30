@doc raw"""
    compute_bergman_time(T::Tracker)

Compute the Bergman time of the mixed cell `σ` with tracker `T`. This is the smallest time t at which `w + t * u` induces a coarser chain of flats than the current chain of flats, where `w` and `u` are the tropical intersection point and tropical drift of `T` respectively.

Caches the answer inside `T`.
"""
function compute_bergman_time(T::Tracker, σ::MixedCell; disable_cache::Bool = false)

    chainOfFlats = chain_of_flats(σ)
    w, u = tropical_intersection_point_and_drift(T, σ)
    inequalities, equalities = cone(chainOfFlats)

    if AbstractAlgebra.get_assertion_level(:TropicalHomotopiesBergman)>0
        # check that we are inside the cone
        if length(equalities) > 0
            @assert all([sum(w .* v) == 0 for v in equalities]) "The intersection point is not in the cone (equality violated)"
            @assert all([sum(w .* v) <= 0 for v in inequalities]) "The intersection point is not in the cone (inequality violated) intersection point: $(w) chain of flats: $(chainOfFlats)"
        else
            @assert all([sum(w .* v) <= 0 for v in inequalities]) "The intersection point is not in the cone"
        end
    end

    if AbstractAlgebra.get_assertion_level(:TropicalHomotopiesBergman)>1
        @assert u in cone_from_equations(linear_equation_matrix(linear_span(C))) "The drift is not in the cone"
    end

    timesOfIntersection = [sum(v.*u) != 0 ? -sum(v.*w) / sum(v.*u) : Nemo.PosInf() for v in inequalities]
    # delete all times that are less than 0
    timesOfIntersection = [t for t in timesOfIntersection if t > 0]
    if timesOfIntersection == []
        return Nemo.PosInf()
    end

    # check that we are still inside the cone
    t = minimum(timesOfIntersection)
    if !isinf(t)
        if length(equalities) > 0
            @assert all([sum((w + t*u) .* v) == 0 for v in equalities]) "The intersection point is not in the cone (equality violated)"
        end
        for v in inequalities
            if sum((w + t*u) .* v) > 0
                @assert false "Inequality $(v) is violated at $(w + t*u) corresponding to t = $(t). The time of intersection for this facet is equal to $(sum(v.*u) != 0 ? -sum(v.*w) / sum(v.*u) : Nemo.PosInf())"
            end
        end
        @assert all([sum((w + t*u) .* v) <= 0 for v in inequalities]) "The intersection point is not in the cone (inequality violated) intersection point: $(w) chain of flats: $(chainOfFlats)"
    end

    bergmanTime = minimum(timesOfIntersection)

    # cache the answer
    if !disable_cache
        update_bergman_time!(T, σ, bergmanTime)
    end

    return bergmanTime
end

@doc raw"""
    bergman_flip(T::Tracker)

Compute the Bergman flip of the mixed cell `σ` with tracker `T`.
"""
function bergman_flip(T::Tracker, σ::MixedCell, tBergman::Height)

    C = chain_of_flats(σ)
    Δ = ambient_support(T)
    w, u = tropical_intersection_point_and_drift(T, σ)

    unrefinedChain = chain_of_flats(matroid(C), w + tBergman * u)

    if length(unrefinedChain) + 1 != length(C)
        return MixedCell[] # need to perturb
    end

    refinedChains = maximal_refinements(unrefinedChain)
    # make sure not to include the original chain
    # refinedChains = [chain for chain in refinedChains if chain != C]

    rows = Vector{QQFieldElem}[]
    for S in supports(active_support(σ))
        pts = points(S)
        p1 = first(pts)
        for p in pts
            if !isequal(p1, p)
                push!(rows, convert(Vector{QQFieldElem}, p1 - p))
            end
        end
    end
    M = Oscar.matrix(QQ, rows)

    jensenTrail = convex_hull([w + tBergman * u], [u], transpose(kernel(M, side=:right)))

    allowedChains = ChainOfFlats[]
    for chain in refinedChains
        # add columns for each indicator vector of chain
        cols = Vector{QQFieldElem}[]
        push!(cols, indicator_vector.(full_flats(chain))...)
        # remove all zero vector from cols
        cols = [col for col in cols if col != zeros(QQ, length(col))]
        A = transpose(Oscar.matrix(QQ, cols))

        # create the matrix formed by the supports
        if Oscar.rank(A) != Oscar.rank(M*A)
            continue
        end

        chainOfFlatsCone = polyhedron(positive_hull(transpose(A), ones_matrix(QQ, 1,nrows(A))))
        if Oscar.dim(intersect(jensenTrail, chainOfFlatsCone)) != 1
            continue
        end

        push!(allowedChains, chain)

    end

    if AbstractAlgebra.get_assertion_level(:TropicalHomotopiesBergman)>0
        # check that the mixed cell data is valid
        newMixedCells = mixed_cell.(Ref(active_support(σ)), allowedChains)
        for σ in newMixedCells
            # check that the matrix coming from σ is invertible
            @assert is_transverse(σ) "$(σ) is not transverse"
            @assert are_support_heights_finite(T, σ) "$(σ) has invalid mixed height data"
        end
        @assert length(newMixedCells) > 0 "No new mixed cells during a Bergman flip"
    end

    return mixed_cell.(Ref(active_support(σ)), allowedChains)

end

@doc raw"""
    bergman_move!(T::Tracker)

Perform a Bergman move on the tracker `T`. This updates the mixed heights and the mixed cells tracked by `T` to the point where a maximal cone of the fine structure is breached, assuming no maximal tropical polyhedra corresponding to the hypersurface dual cells are breached.

We can perform a Bergman move when `bergman_time(T)` is less than `jensen_time(T)`.
"""
function bergman_move!(T::Tracker)

    # work out which mixed cell(s) have minimal Bergman time
    smallestTBergman = minimum([bergman_time(T, σ) for σ in mixed_cells(T)])
    changingMixedCells = [σ for σ in mixed_cells(T) if bergman_time(T, σ) == smallestTBergman]

    newMixedCells = MixedCell[]
    for σ in changingMixedCells
        push!(newMixedCells, bergman_flip(T, σ)...)
        remove_mixed_cell!(T, σ)
    end

    # move the tracker to the heights achieved at the bergman time
    add_heights!(T, smallestTBergman*direction(T))

    for σ in newMixedCells
        merge_mixed_cell!(T, σ)
    end

end

# orthogonal projection away from subspace spanned by columns of kernelMatrix
function projection_matrix(kernelMatrix)
    _, invKernelMatrix = Oscar.is_invertible_with_inverse(transpose(kernelMatrix) * kernelMatrix)
    return kernelMatrix * invKernelMatrix * transpose(kernelMatrix)
end

@doc raw"""
    oblique_projection_matrix(A::Matrix, B::Matrix)

Compute the oblique projection matrix from the subspace spanned by the columns of `A` onto the subspace spanned by the columns of `B`.
"""
function oblique_projection_matrix(A, B)
    if AbstractAlgebra.get_assertion_level(:TropicalHomotopiesBergman)>0
        @assert Oscar.rank(B) == ncols(B) "Matrix B is not full rank"
        @assert Oscar.rank(A) == ncols(A) "Matrix A is not full rank"
    end
    _, invKernelMatrix = Oscar.is_invertible_with_inverse(transpose(B) * A)
    return A*(invKernelMatrix)*transpose(B)
end
