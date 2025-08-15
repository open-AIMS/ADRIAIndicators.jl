"""
    _reef_condition_index!(rc::AbstractArray, evenness::AbstractArray, sv::AbstractArray, juves::AbstractArray, out_rci::AbstractArray; threshold=2)::Nothing

Calculate the Reef Condition Index (RCI).

# Arguments
- `rc` : Relative coral cover.
- `evenness` : Coral evenness.
- `sv` : Relative shelter volume.
- `out_rci` : Output array buffer for the RCI.
"""
function _reef_condition_index!(
    rc::AbstractArray{T,N},
    cd::AbstractArray{T,N},
    sv::AbstractArray{T,N},
    out_rci::AbstractArray{T,N}
)::Nothing where {T<:AbstractFloat,N}
    out_rci .= clamp.((rc .+ cd .+ sv) ./ 3.0, 0.0, 1.0)
    return nothing
end

"""
    reef_condition_index(rc::AbstractArray, evenness::AbstractArray, sv::AbstractArray, juves::AbstractArray; threshold=2)::AbstractArray

Calculate the Reef Condition Index (RCI).

# Arguments
- `rc` : Relative coral cover.
- `evenness` : Coral evenness.
- `sv` : Relative shelter volume.

# Returns
A 2D array of the Reef Condition Index.
"""
function reef_condition_index(
    rc::AbstractArray,
    evenness::AbstractArray,
    sv::AbstractArray
)::AbstractArray
    if (size(rc) != size(evenness)) || (size(rc) != size(sv))
        throw(DimensionMismatch("All input metric arrays must have the same dimensions."))
    end

    out_rci = zeros(Float64, size(rc))
    _reef_condition_index!(rc, evenness, sv, out_rci)

    return out_rci
end

"""
    _reef_tourism_index!(rc::AbstractArray, evenness::AbstractArray, sv::AbstractArray, juves::AbstractArray, intcp_u::Real, out_rti::AbstractArray)::Nothing

Calculate the Reef Tourism Index (RTI) for a single scenario.

# Arguments
- `rc` : Relative coral cover.
- `evenness` : Coral evenness.
- `sv` : Relative shelter volume.
- `juves` : Abundance of coral juveniles.
- `intcp_u` : A scenario-specific intercept.
- `out_rti` : Output array buffer for the RTI.
"""
function _reef_tourism_index!(
    rc::AbstractArray,
    evenness::AbstractArray,
    sv::AbstractArray,
    juves::AbstractArray,
    intcp_u::Real,
    out_rti::AbstractArray
)::Nothing
    intcp = 0.47947 + intcp_u

    out_rti .= intcp .+ (0.12764 .* rc) .+ (0.31946 .* evenness) .+ (0.11676 .* sv) .+ (-0.0036065 .* juves)
    out_rti .= round.(clamp.(out_rti, 0.1, 0.9), digits=2)

    return nothing
end

"""
    reef_tourism_index(rc::AbstractArray, evenness::AbstractArray, sv::AbstractArray, juves::AbstractArray, intcp_u::Real)::AbstractArray

Calculate the Reef Tourism Index (RTI) for a single scenario.

# Arguments
- `rc` : Relative coral cover.
- `evenness` : Coral evenness.
- `sv` : Relative shelter volume.
- `juves` : Abundance of coral juveniles.
- `intcp_u` : A scenario-specific intercept.

# Returns
A 2D array of the Reef Tourism Index.
"""
function reef_tourism_index(
    rc::AbstractArray,
    evenness::AbstractArray,
    sv::AbstractArray,
    juves::AbstractArray,
    intcp_u::Real
)::AbstractArray
    if (size(rc) != size(evenness)) || (size(rc) != size(sv)) || (size(rc) != size(juves))
        throw(DimensionMismatch("All input metric arrays must have the same dimensions."))
    end

    out_rti = zeros(Float64, size(rc))
    _reef_tourism_index!(rc, evenness, sv, juves, intcp_u, out_rti)

    return out_rti
end

"""
    _reef_fish_index!(rc::AbstractArray, intcp_u1::Real, intcp_u2::Real, out_rfi::AbstractArray)::Nothing

Calculate the Reef Fish Index (RFI) for a single scenario.

# Arguments
- `rc` : Relative coral cover.
- `intcp_u1` : Scenario-specific intercept for the first linear model.
- `intcp_u2` : Scenario-specific intercept for the second linear model.
- `out_rfi` : Output array buffer for the RFI.
"""
function _reef_fish_index!(
    rc::AbstractArray,
    intcp_u1::Real,
    intcp_u2::Real,
    out_rfi::AbstractArray
)::Nothing
    intcp1 = 1.232 + intcp_u1
    intcp2 = -1623.6 + intcp_u2
    slope1 = 0.007476
    slope2 = 1883.3

    out_rfi .= 0.01 .* (intcp2 .+ slope2 .* (intcp1 .+ slope1 .* (rc .* 100.0)))
    out_rfi .= round.(out_rfi, digits=2)

    return nothing
end

"""
    reef_fish_index(rc::AbstractArray, intcp_u1::Real, intcp_u2::Real)::AbstractArray

Calculate the Reef Fish Index (RFI) for a single scenario.

# Arguments
- `rc` : Relative coral cover.
- `intcp_u1` : Scenario-specific intercept for the first linear model.
- `intcp_u2` : Scenario-specific intercept for the second linear model.

# Returns
A 2D array of the Reef Fish Index.
"""
function reef_fish_index(
    rc::AbstractArray,
    intcp_u1::Real,
    intcp_u2::Real
)::AbstractArray
    out_rfi = zeros(Float64, size(rc))
    _reef_fish_index!(rc, intcp_u1, intcp_u2, out_rfi)
    return out_rfi
end
