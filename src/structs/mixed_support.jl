using Oscar
export MixedSupport, vector_of_points, cayley_embedding, mixed_subdivision, supports

@doc raw"""
    struct MixedSupport

A mixed support is a collection of supports. This can represent either ambient support or a mixed cell candidate.
"""
mutable struct MixedSupport

    supports::Tuple{Vararg{Support}} # collection of supports

end

function mixed_support(S::Tuple{Vararg{Support}})::MixedSupport
    return MixedSupport(S)
end

function mixed_support()::MixedSupport
    return mixed_support(tuple())
end

function supports(Δ::MixedSupport)::Tuple{Vararg{Support}}
    return Δ.supports
end

function Base.length(Δ::MixedSupport)
    return length(supports(Δ))
end

function Base.show(io::IO, Δ::MixedSupport)

    print(io, "Mixed support involving $(length(supports(Δ))) point supports")
end

@doc raw"""
    vector_of_points(Δ::MixedSupport)

Convenience function to return the points of a mixed support as a partitioned vector.
"""
function vector_of_points(Δ::MixedSupport)
    return [[p for p in points(s)] for s in supports(Δ)]
end

function points(Δ::MixedSupport)
    return collect(Iterators.flatten(vector_of_points(Δ)))
end

@doc raw"""
    mixed_subdivision(Δ::MixedSupport)

Computes the mixed subdivision corresponding to the mixed support Δ.
"""
function mixed_subdivision(Δ::MixedSupport)
    minkowskiSum = minkowski_sum(points.(supports(Δ))...)
    mixedHeights = minkowski_sum(heights.(supports(Δ))...)

    subdivision = subdivision_of_points(convert.(Vector{Int}, minkowskiSum), convert.(QQFieldElem, mixedHeights))
    return subdivision
end

@doc raw"""
    is_subset(S::MixedSupport, T::MixedSupport)

Returns `true if `S` is a subset of `T`, and `false` otherwise.
"""
function is_subset(S::MixedSupport, T::MixedSupport)

    if length(S) != length(T)
        return false
    end

    for (s, t) in zip(supports(S), supports(T))
        if !is_subset(s, t)
            return false
        end
    end

    return true
end

@doc raw"""
    mixed_support(ambientSupport::MixedSupport, oldSupport::Support, newSupport::Support)

Given an ambient support `ambientSupport`, a support `oldSupport` in `ambientSupport`, and a new support `newSupport` that is a superset of `oldSupport`, return the mixed support that is the result of replacing `oldSupport` with `newSupport` in `ambientSupport`.
"""
function mixed_support(ambientSupport::MixedSupport, oldSupport::Support, newSupport::Support)

    @assert is_subset(oldSupport, newSupport) "The new support must be a superset of the old support"

    newSupports = [s == oldSupport ? newSupport : s for s in supports(ambientSupport)]
    return mixed_support(tuple(newSupports...))
end

function Base.getindex(Δ::MixedSupport, p::Point)
    index = findfirst(s -> in(p, s), supports(Δ))
    if isnothing(index)
        return Nemo.PosInf()
    end
    memberSupport = supports(Δ)[findfirst(s -> in(p, s), supports(Δ))]
    @assert !isnothing(memberSupport) "The point $p is not in the mixed support"
    return memberSupport[p]
end

function Base.:-(a::MixedSupport, b::MixedSupport)::MixedSupport
    return mixed_support(tuple([s - t for (s, t) in zip(supports(a), supports(b))]...))
end

function Base.:*(a::Height, Δ::MixedSupport)::MixedSupport
    return mixed_support(tuple([a * s for s in supports(Δ)]...))
end

function combine(Δ::MixedSupport, Ε::MixedSupport)
    return mixed_support(tuple(supports(Δ)..., supports(Ε)...))
end