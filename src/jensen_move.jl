export jensen_time

@doc raw"""
    jensen_time(T::Tracker, σ::MixedCell)

Return the time at which the mixed cell candidate `σ` breaches its mixed cell cone.
"""
function jensen_time(T::Tracker, σ::MixedCell)::Union{QQFieldElem,PosInf}
    
    hypersurfaceDuals = transform_linear_support(chain_of_flats(σ))

    @assert is_subset(active_support(σ), ambient_support(T)) "The active support of the mixed cell is not a subset of the ambient support."
    δ = combine(active_support(σ), hypersurfaceDuals)
    Δ = combine(ambient_support(T), hypersurfaceDuals)

    C = mixed_cell_cone(δ, Δ)

    @assert Δ in C "The mixed cell being tracked is not in the mixed cell cone."

    v = direction(T)

    timesOfIntersection = [dot(v, κ) != 0 ? -dot(κ, Δ) / dot(v, κ) : Nemo.PosInf() for κ in facets(C)]
    # delete all times that are less than 0
    timesOfIntersection = [t for t in timesOfIntersection if t > 0]
    if timesOfIntersection == []
        return Nemo.PosInf()
    end

    return minimum(timesOfIntersection)
    
end


function jensen_flip(T::Tracker, σ::MixedCell, tJensen::Height)

    hypersurfaceDuals = transform_linear_support(chain_of_flats(σ))

    δ = combine(active_support(σ), hypersurfaceDuals)
    Δ = combine(ambient_support(T), hypersurfaceDuals)

    v = direction(T)

    C = mixed_cell_cone(δ, Δ)

    @assert sum([dot(v, κ) * tJensen == -dot(Δ, κ) for κ in facets(C)]) == 1 "The mixed cell being tracked does not breach its mixed cell cone in exactly one facet."

    κ = facets(C)[findfirst(κ -> dot(v, κ) * tJensen == -dot(Δ, κ), facets(C))]
    @debug "The facet inequality broken is given by $(κ)"
    p = extra_point(κ)
    changingSupport = supports(Δ)[findfirst(p in points(S) for S in supports(Δ))]
    # work out which mixed cell support this corresponds to (which active support)
    changingDualCell = supports(σ)[findfirst(is_subset(S, changingSupport) for S in supports(σ))]

    @assert length(changingDualCell) == 2 "The changing dual cell is not minimal."
    newMixedCells = MixedCell[]

    for q in points(changingDualCell)
        if κ[q] > 0
            @debug "Adding mixed cell with $q in support and $p not in support"
            # return mixed cell with p in support and q not in support
            push!(newMixedCells, swap(σ, q, p))
        end
    end

    for σ in newMixedCells
        # check that the matrix coming from σ is invertible
        @assert is_transverse(σ) "$(σ) is not transverse"
        @assert are_support_heights_finite(T, σ) "$(σ) has invalid mixed height data"

    end

    return newMixedCells
end

@doc raw"""
    jensen_move!(T::Tracker)

Perform a Jensen move on the tracker `T`.
"""
function jensen_move!(T::Tracker)
    
    # work out which mixed cell(s) have minimal Bergman time
    smallestTJensen = minimum([jensen_time(T, σ) for σ in mixed_cells(T)])
    changingMixedCells = [σ for σ in mixed_cells(T) if jensen_time(T, σ) == smallestTJensen]

    newMixedCells = MixedCell[]
    for σ in changingMixedCells
        push!(newMixedCells, jensen_flip(T, σ)...)
        remove_mixed_cell!(T, σ)
    end

    # move the tracker to the heights achieved at the bergman time
    add_heights!(T, smallestTJensen*direction(T))

    for σ in newMixedCells
        merge_mixed_cell!(T, σ)
    end

end