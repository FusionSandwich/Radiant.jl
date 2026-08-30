"""
    Energy_Group_Condensation_Map

Exact nested mapping from fine ascending energy groups to coarser ascending energy groups. Every
coarse boundary must coincide with one fine boundary and the full energy interval must be retained.
"""
struct Energy_Group_Condensation_Map
    fine_edges_eV::Vector{Float64}
    coarse_edges_eV::Vector{Float64}
    fine_to_coarse::Vector{Int64}
    coarse_members::Vector{Vector{Int64}}

    function Energy_Group_Condensation_Map(
        fine_edges_eV::AbstractVector{<:Real},
        coarse_edges_eV::AbstractVector{<:Real};
        edge_rtol::Real=1.0e-12,
        edge_atol_eV::Real=1.0e-10,
    )
        fine = Float64.(fine_edges_eV)
        coarse = Float64.(coarse_edges_eV)
        length(fine) >= 2 && all(isfinite,fine) && all(fine .>= 0.0) &&
            all(diff(fine) .> 0.0) || error(
            "Fine energy boundaries must be finite, nonnegative, and strictly increasing.",
        )
        length(coarse) >= 2 && all(isfinite,coarse) && all(coarse .>= 0.0) &&
            all(diff(coarse) .> 0.0) || error(
            "Coarse energy boundaries must be finite, nonnegative, and strictly increasing.",
        )
        isapprox(coarse[1],fine[1];rtol=edge_rtol,atol=edge_atol_eV) &&
            isapprox(coarse[end],fine[end];rtol=edge_rtol,atol=edge_atol_eV) || error(
            "Coarse groups must preserve the complete fine energy interval.",
        )
        fine_boundary_indices = Int64[]
        for boundary in coarse
            matches = findall(value -> isapprox(
                value,boundary;rtol=edge_rtol,atol=edge_atol_eV,
            ),fine)
            length(matches) == 1 || error(
                "Every coarse boundary must match exactly one fine boundary.",
            )
            push!(fine_boundary_indices,Int64(first(matches)))
        end
        all(diff(fine_boundary_indices) .> 0) || error(
            "Coarse boundary mapping is not strictly increasing.",
        )
        fine_to_coarse = zeros(Int64,length(fine)-1)
        members = Vector{Vector{Int64}}(undef,length(coarse)-1)
        for coarse_group in 1:length(coarse)-1
            first_fine = fine_boundary_indices[coarse_group]
            last_fine = fine_boundary_indices[coarse_group+1]-1
            first_fine <= last_fine || error("A coarse group contains no fine groups.")
            members[coarse_group] = collect(Int64(first_fine):Int64(last_fine))
            fine_to_coarse[first_fine:last_fine] .= coarse_group
        end
        all(fine_to_coarse .> 0) || error("At least one fine group was not mapped.")
        return new(fine,coarse,fine_to_coarse,members)
    end
end

fine_group_count(mapping::Energy_Group_Condensation_Map) = length(mapping.fine_to_coarse)
coarse_group_count(mapping::Energy_Group_Condensation_Map) = length(mapping.coarse_members)

function condense_group_integrals(
    mapping::Energy_Group_Condensation_Map,
    fine_integrals::AbstractVector{<:Real},
)
    fine = Float64.(fine_integrals)
    length(fine) == fine_group_count(mapping) || error(
        "Fine group-integral vector does not match the condensation map.",
    )
    all(isfinite,fine) || error("Fine group integrals must be finite.")
    coarse = zeros(Float64,coarse_group_count(mapping))
    for fine_group in eachindex(fine)
        coarse[mapping.fine_to_coarse[fine_group]] += fine[fine_group]
    end
    return coarse
end

"""
Return a coarse-by-fine linear weighting matrix. For each coarse group, the fine reference fluxes
inside that group are normalized to one. A zero-flux coarse group is rejected; supplying a response
coefficient for an unsampled energy interval would otherwise be arbitrary.
"""
function response_condensation_weight_matrix(
    mapping::Energy_Group_Condensation_Map,
    fine_reference_flux::AbstractVector{<:Real},
)
    flux = Float64.(fine_reference_flux)
    length(flux) == fine_group_count(mapping) || error(
        "Reference flux does not match the fine energy groups.",
    )
    all(value -> isfinite(value) && value >= 0.0,flux) || error(
        "Reference group flux must be finite and nonnegative.",
    )
    weights = zeros(Float64,coarse_group_count(mapping),fine_group_count(mapping))
    for coarse_group in 1:coarse_group_count(mapping)
        members = mapping.coarse_members[coarse_group]
        denominator = sum(flux[members])
        denominator > 0.0 || error(
            "Coarse group $(coarse_group) has zero reference flux; a protected-response " *
            "coefficient cannot be generated without an explicit alternate spectrum.",
        )
        weights[coarse_group,members] .= flux[members]./denominator
    end
    return weights
end

struct Response_Condensation_Receipt
    response_id::String
    fine_response_integral::Float64
    coarse_response_integral::Float64
    absolute_residual::Float64
    relative_residual::Float64
    reference_flux_hash::String
    group_map_hash::String
    physical_reference::Bool
    metadata::Dict{String,String}
end

function _numeric_vector_hash(values::AbstractVector{<:Real})
    io = IOBuffer()
    for value in Float64.(values)
        Base.write(io,repr(value))
        Base.write(io,UInt8('\n'))
    end
    return bytes2hex(SHA.sha256(take!(io)))
end

function group_condensation_map_hash(mapping::Energy_Group_Condensation_Map)
    io = IOBuffer()
    Base.write(io,"radiant.hts.energy_group_condensation_map/v1\n")
    for edge in mapping.fine_edges_eV
        Base.write(io,repr(edge)); Base.write(io,UInt8('\n'))
    end
    Base.write(io,"coarse\n")
    for edge in mapping.coarse_edges_eV
        Base.write(io,repr(edge)); Base.write(io,UInt8('\n'))
    end
    return bytes2hex(SHA.sha256(take!(io)))
end

"""
    condense_response_coefficients(mapping, flux, response; response_id, ...)

Flux-weight one fine-group response coefficient into each coarse group. The resulting coefficient
preserves `sum(flux .* response)` exactly for the declared reference spectrum, up to floating-point
roundoff. It does not claim preservation for a different spectrum.
"""
function condense_response_coefficients(
    mapping::Energy_Group_Condensation_Map,
    fine_reference_flux::AbstractVector{<:Real},
    fine_response_coefficients::AbstractVector{<:Real};
    response_id::AbstractString,
    physical_reference::Bool=false,
    closure_rtol::Real=1.0e-12,
    closure_atol::Real=1.0e-12,
    metadata::AbstractDict=Dict{String,String}(),
)
    isempty(response_id) && error("Condensed response identifier cannot be empty.")
    flux = Float64.(fine_reference_flux)
    response = Float64.(fine_response_coefficients)
    length(response) == fine_group_count(mapping) || error(
        "Fine response coefficients do not match the condensation map.",
    )
    all(isfinite,response) || error("Fine response coefficients must be finite.")
    weights = response_condensation_weight_matrix(mapping,flux)
    coarse_response = weights*response
    coarse_flux = condense_group_integrals(mapping,flux)
    fine_integral = dot(flux,response)
    coarse_integral = dot(coarse_flux,coarse_response)
    residual = coarse_integral-fine_integral
    relative = abs(residual)/max(abs(fine_integral),abs(coarse_integral),eps(Float64))
    isapprox(coarse_integral,fine_integral;rtol=closure_rtol,atol=closure_atol) || error(
        "Response-preserving condensation failed closure for $(response_id).",
    )
    metadata_string = Dict{String,String}(
        "weighting" => "reference-flux",
        "spectrum_specific" => "true",
        "fine_group_count" => string(fine_group_count(mapping)),
        "coarse_group_count" => string(coarse_group_count(mapping)),
    )
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    receipt = Response_Condensation_Receipt(
        String(response_id),fine_integral,coarse_integral,abs(residual),relative,
        _numeric_vector_hash(flux),group_condensation_map_hash(mapping),physical_reference,
        metadata_string,
    )
    return coarse_response,receipt
end

"""Propagate a full fine-group covariance matrix through the spectrum-dependent linear weights."""
function condense_response_covariance(
    mapping::Energy_Group_Condensation_Map,
    fine_reference_flux::AbstractVector{<:Real},
    fine_covariance::AbstractMatrix{<:Real};
    symmetry_tolerance::Real=1.0e-10,
)
    covariance = Float64.(fine_covariance)
    expected = fine_group_count(mapping)
    size(covariance) == (expected,expected) || error(
        "Fine response covariance must have shape (Nfine,Nfine).",
    )
    all(isfinite,covariance) || error("Fine response covariance must be finite.")
    isapprox(covariance,covariance';rtol=symmetry_tolerance,atol=symmetry_tolerance) || error(
        "Fine response covariance must be symmetric.",
    )
    minimum(diag(covariance)) >= -Float64(symmetry_tolerance) || error(
        "Fine response covariance has a negative diagonal variance.",
    )
    weights = response_condensation_weight_matrix(mapping,fine_reference_flux)
    condensed = weights*covariance*weights'
    return 0.5.*(condensed+condensed')
end

struct Transfer_Matrix_Condensation_Receipt
    response_id::String
    fine_outgoing_source::Vector{Float64}
    coarse_outgoing_source::Vector{Float64}
    aggregated_fine_outgoing_source::Vector{Float64}
    maximum_absolute_residual::Float64
    maximum_relative_residual::Float64
    incoming_flux_hash::String
    incoming_group_map_hash::String
    outgoing_group_map_hash::String
    physical_reference::Bool
    metadata::Dict{String,String}
end

"""
    condense_transfer_matrix(in_map, out_map, fine_flux, fine_transfer; ...)

Condense a transfer/production matrix stored as `(incoming_group,outgoing_group)`. The incoming
reference flux weights each incoming coarse row, while outgoing fine groups are summed into their
coarse group. This exactly preserves the aggregated outgoing source for the reference spectrum.
It is suitable for neutron-to-photon production, scattering yields, and response matrices when the
matrix orientation and units are explicitly declared.
"""
function condense_transfer_matrix(
    incoming_mapping::Energy_Group_Condensation_Map,
    outgoing_mapping::Energy_Group_Condensation_Map,
    fine_incoming_flux::AbstractVector{<:Real},
    fine_transfer::AbstractMatrix{<:Real};
    response_id::AbstractString,
    physical_reference::Bool=false,
    closure_rtol::Real=1.0e-12,
    closure_atol::Real=1.0e-12,
    metadata::AbstractDict=Dict{String,String}(),
)
    isempty(response_id) && error("Condensed transfer response identifier cannot be empty.")
    flux = Float64.(fine_incoming_flux)
    matrix = Float64.(fine_transfer)
    size(matrix) == (
        fine_group_count(incoming_mapping),fine_group_count(outgoing_mapping),
    ) || error(
        "Fine transfer matrix must have shape (fine incoming groups,fine outgoing groups).",
    )
    all(isfinite,matrix) || error("Fine transfer matrix must be finite.")
    weights = response_condensation_weight_matrix(incoming_mapping,flux)
    coarse_matrix = zeros(
        Float64,coarse_group_count(incoming_mapping),coarse_group_count(outgoing_mapping),
    )
    for incoming_coarse in 1:coarse_group_count(incoming_mapping)
        for outgoing_coarse in 1:coarse_group_count(outgoing_mapping)
            outgoing_members = outgoing_mapping.coarse_members[outgoing_coarse]
            for incoming_fine in incoming_mapping.coarse_members[incoming_coarse]
                coarse_matrix[incoming_coarse,outgoing_coarse] +=
                    weights[incoming_coarse,incoming_fine]*
                    sum(matrix[incoming_fine,outgoing_members])
            end
        end
    end
    fine_outgoing = vec(transpose(flux)*matrix)
    aggregated_fine = condense_group_integrals(outgoing_mapping,fine_outgoing)
    coarse_flux = condense_group_integrals(incoming_mapping,flux)
    coarse_outgoing = vec(transpose(coarse_flux)*coarse_matrix)
    residual = coarse_outgoing-aggregated_fine
    scale = max.(max.(abs.(coarse_outgoing),abs.(aggregated_fine)),eps(Float64))
    relative = abs.(residual)./scale
    isapprox(coarse_outgoing,aggregated_fine;rtol=closure_rtol,atol=closure_atol) || error(
        "Transfer-matrix condensation failed outgoing-source closure for $(response_id).",
    )
    metadata_string = Dict{String,String}(
        "matrix_orientation" => "incoming-group-by-outgoing-group",
        "incoming_weighting" => "reference-flux",
        "outgoing_operation" => "group-sum",
        "spectrum_specific" => "true",
    )
    for (key,value) in metadata
        metadata_string[string(key)] = string(value)
    end
    receipt = Transfer_Matrix_Condensation_Receipt(
        String(response_id),fine_outgoing,coarse_outgoing,aggregated_fine,
        maximum(abs.(residual)),maximum(relative),_numeric_vector_hash(flux),
        group_condensation_map_hash(incoming_mapping),
        group_condensation_map_hash(outgoing_mapping),physical_reference,metadata_string,
    )
    return coarse_matrix,receipt
end

function response_condensation_receipt(receipt::Response_Condensation_Receipt)
    return Dict{String,Any}(
        "schema" => "radiant.hts.response_condensation_receipt/v1",
        "response_id" => receipt.response_id,
        "fine_response_integral" => receipt.fine_response_integral,
        "coarse_response_integral" => receipt.coarse_response_integral,
        "absolute_residual" => receipt.absolute_residual,
        "relative_residual" => receipt.relative_residual,
        "reference_flux_hash" => receipt.reference_flux_hash,
        "group_map_hash" => receipt.group_map_hash,
        "physical_reference" => receipt.physical_reference,
        "metadata" => copy(receipt.metadata),
    )
end

function transfer_matrix_condensation_receipt(
    receipt::Transfer_Matrix_Condensation_Receipt,
)
    return Dict{String,Any}(
        "schema" => "radiant.hts.transfer_matrix_condensation_receipt/v1",
        "response_id" => receipt.response_id,
        "fine_outgoing_source" => receipt.fine_outgoing_source,
        "coarse_outgoing_source" => receipt.coarse_outgoing_source,
        "aggregated_fine_outgoing_source" => receipt.aggregated_fine_outgoing_source,
        "maximum_absolute_residual" => receipt.maximum_absolute_residual,
        "maximum_relative_residual" => receipt.maximum_relative_residual,
        "incoming_flux_hash" => receipt.incoming_flux_hash,
        "incoming_group_map_hash" => receipt.incoming_group_map_hash,
        "outgoing_group_map_hash" => receipt.outgoing_group_map_hash,
        "physical_reference" => receipt.physical_reference,
        "metadata" => copy(receipt.metadata),
    )
end
