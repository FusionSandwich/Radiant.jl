struct Coupling_Solver_Output
    boundary_source::Boundary_Angular_Current_Source
    protected_responses::Vector{Float64}
    particle_balance_residual::Float64
    energy_balance_residual::Float64
    clipping_used::Bool
    negative_source_repair_used::Bool
    solver_id::String
    result_artifact_hash::String
    provenance::Dict{String,String}

    function Coupling_Solver_Output(
        boundary_source::Boundary_Angular_Current_Source;
        protected_responses::AbstractVector{<:Real},
        particle_balance_residual::Real,
        energy_balance_residual::Real,
        clipping_used::Bool=false,
        negative_source_repair_used::Bool=false,
        solver_id::AbstractString,
        result_artifact_hash::AbstractString,
        provenance::AbstractDict=Dict{String,String}(),
    )
        responses = Float64.(protected_responses)
        isempty(responses) && error("Coupling solver output requires protected responses.")
        all(isfinite,responses) || error("Coupling protected responses must be finite.")
        particle_residual = Float64(particle_balance_residual)
        energy_residual = Float64(energy_balance_residual)
        all(isfinite,(particle_residual,energy_residual)) || error(
            "Coupling balance residuals must be finite.",
        )
        isempty(solver_id) && error("Coupling solver identifier cannot be empty.")
        isempty(result_artifact_hash) && error("Coupling result hash cannot be empty.")
        provenance_string = Dict{String,String}()
        for (key,value) in provenance
            provenance_string[string(key)] = string(value)
        end
        return new(
            boundary_source,responses,particle_residual,energy_residual,clipping_used,
            negative_source_repair_used,String(solver_id),String(result_artifact_hash),
            provenance_string,
        )
    end
end

struct Closed_Coupling_Settings
    maximum_iterations::Int64
    minimum_iterations::Int64
    relaxation::Float64
    current_relative_tolerance::Float64
    energy_current_relative_tolerance::Float64
    response_relative_tolerance::Float64
    particle_balance_tolerance::Float64
    energy_balance_tolerance::Float64

    function Closed_Coupling_Settings(;
        maximum_iterations::Integer=50,
        minimum_iterations::Integer=2,
        relaxation::Real=0.5,
        current_relative_tolerance::Real=1.0e-6,
        energy_current_relative_tolerance::Real=1.0e-6,
        response_relative_tolerance::Real=1.0e-5,
        particle_balance_tolerance::Real=1.0e-8,
        energy_balance_tolerance::Real=1.0e-8,
    )
        maximum_iterations >= 1 || error("Maximum coupling iterations must be positive.")
        minimum_iterations >= 1 && minimum_iterations <= maximum_iterations || error(
            "Minimum coupling iterations must lie inside the iteration range.",
        )
        values = Float64[
            relaxation,current_relative_tolerance,energy_current_relative_tolerance,
            response_relative_tolerance,particle_balance_tolerance,energy_balance_tolerance,
        ]
        all(isfinite,values) || error("Coupling settings must be finite.")
        0.0 < values[1] <= 1.0 || error("Coupling relaxation must lie in (0,1].")
        all(values[2:end] .>= 0.0) || error("Coupling tolerances must be nonnegative.")
        return new(Int64(maximum_iterations),Int64(minimum_iterations),values...)
    end
end

struct Coupling_Iteration_Record
    iteration::Int64
    forward_current_change::Float64
    forward_energy_current_change::Float64
    return_current_change::Float64
    return_energy_current_change::Float64
    protected_response_change::Float64
    opensn_particle_balance_residual::Float64
    opensn_energy_balance_residual::Float64
    radiant_particle_balance_residual::Float64
    radiant_energy_balance_residual::Float64
    converged::Bool
end

struct Closed_Coupling_Result
    forward_source::Boundary_Angular_Current_Source
    return_source::Boundary_Angular_Current_Source
    opensn_output::Coupling_Solver_Output
    radiant_output::Coupling_Solver_Output
    iterations::Vector{Coupling_Iteration_Record}
    converged::Bool
    ownership_map_hash::String
    forward_source_hash::String
    return_source_hash::String
    provenance::Dict{String,String}
end

function _boundary_sources_compatible(
    first::Boundary_Angular_Current_Source,
    second::Boundary_Angular_Current_Source;
    tolerance::Real=1.0e-12,
)
    get_tag(first.particle) == get_tag(second.particle) || return false
    first.patch_ids == second.patch_ids || return false
    first.energy_edges_eV == second.energy_edges_eV || return false
    first.quadrature_weights == second.quadrature_weights || return false
    size(first.angular_flux) == size(second.angular_flux) || return false
    isapprox(first.centroids_cm,second.centroids_cm;rtol=tolerance,atol=tolerance) || return false
    isapprox(first.areas_cm2,second.areas_cm2;rtol=tolerance,atol=tolerance) || return false
    isapprox(first.normals,second.normals;rtol=tolerance,atol=tolerance) || return false
    isapprox(first.directions,second.directions;rtol=tolerance,atol=tolerance) || return false
    return true
end

function scale_boundary_source(
    source::Boundary_Angular_Current_Source,
    factor::Real;
    source_hash::AbstractString=source.normalization.source_hash,
)
    value = Float64(factor)
    isfinite(value) && value >= 0.0 || error("Boundary-source scale must be nonnegative.")
    variance = isnothing(source.variance) ? nothing : value^2 .* source.variance
    normalization = Source_Normalization(
        basis=source.normalization.basis,
        source_rate_per_s=source.normalization.source_rate_per_s,
        symmetry_factor=source.normalization.symmetry_factor,
        time_interval_s=source.normalization.time_interval_s,
        time_class=source.normalization.time_class,
        source_hash=source_hash,
        provenance=merge(
            copy(source.normalization.provenance),
            Dict("scaled_by" => string(value)),
        ),
    )
    return Boundary_Angular_Current_Source(
        source.particle,source.patch_ids,source.centroids_cm,source.areas_cm2,source.normals,
        source.tangent_1,source.tangent_2,source.energy_edges_eV,source.directions,
        source.quadrature_weights,value.*source.angular_flux,normalization;
        variance=variance,provenance=merge(copy(source.provenance),Dict(
            "source_operation" => "scale",
            "scale_factor" => string(value),
        )),
    )
end

function blend_boundary_sources(
    previous::Boundary_Angular_Current_Source,
    candidate::Boundary_Angular_Current_Source,
    relaxation::Real;
    source_hash::AbstractString="closed-coupling-blend",
)
    _boundary_sources_compatible(previous,candidate) || error(
        "Boundary sources are incompatible and cannot be relaxed.",
    )
    factor = Float64(relaxation)
    isfinite(factor) && 0.0 < factor <= 1.0 || error(
        "Boundary-source relaxation must lie in (0,1].",
    )
    values = (1.0-factor).*previous.angular_flux+factor.*candidate.angular_flux
    any(values .< 0.0) && error(
        "Relaxed boundary source became negative; clipping is not permitted.",
    )
    variance = if isnothing(previous.variance) && isnothing(candidate.variance)
        nothing
    else
        previous_variance = isnothing(previous.variance) ? zeros(size(values)) : previous.variance
        candidate_variance = isnothing(candidate.variance) ? zeros(size(values)) : candidate.variance
        (1.0-factor)^2 .* previous_variance+factor^2 .* candidate_variance
    end
    normalization = Source_Normalization(
        basis=candidate.normalization.basis,
        source_rate_per_s=candidate.normalization.source_rate_per_s,
        symmetry_factor=candidate.normalization.symmetry_factor,
        time_interval_s=candidate.normalization.time_interval_s,
        time_class=candidate.normalization.time_class,
        source_hash=source_hash,
        provenance=Dict(
            "operation" => "linear-relaxation",
            "relaxation" => string(factor),
            "previous_hash" => previous.normalization.source_hash,
            "candidate_hash" => candidate.normalization.source_hash,
        ),
    )
    return Boundary_Angular_Current_Source(
        candidate.particle,candidate.patch_ids,candidate.centroids_cm,candidate.areas_cm2,
        candidate.normals,candidate.tangent_1,candidate.tangent_2,
        candidate.energy_edges_eV,candidate.directions,candidate.quadrature_weights,values,
        normalization;
        variance=variance,provenance=merge(copy(candidate.provenance),Dict(
            "source_operation" => "closed-coupling-relaxation",
        )),
    )
end

function boundary_energy_current(
    source::Boundary_Angular_Current_Source;
    representative_energy_eV::Union{Nothing,AbstractVector{<:Real}}=nothing,
    physical::Bool=false,
)
    energies = if isnothing(representative_energy_eV)
        0.5 .* (source.energy_edges_eV[1:end-1]+source.energy_edges_eV[2:end])
    else
        Float64.(representative_energy_eV)
    end
    length(energies) == length(source.energy_edges_eV)-1 || error(
        "Representative-energy vector does not match boundary-source groups.",
    )
    all(value -> isfinite(value) && value >= 0.0,energies) || error(
        "Representative energies must be finite and nonnegative.",
    )
    value = sum(get_incoming_current(source).*reshape(energies,1,length(energies)))
    return physical ? apply_normalization(value,source.normalization) : value
end

function _relative_scalar_change(current::Real,previous::Real,atol::Real=1.0e-30)
    return abs(Float64(current)-Float64(previous))/max(abs(Float64(current)),abs(Float64(previous)),Float64(atol))
end

function boundary_current_relative_change(
    current::Boundary_Angular_Current_Source,
    previous::Boundary_Angular_Current_Source,
)
    _boundary_sources_compatible(current,previous) || error(
        "Boundary-current change requires compatible sources.",
    )
    current_values = get_incoming_current(current)
    previous_values = get_incoming_current(previous)
    denominator = max.(abs.(current_values),abs.(previous_values),eps(Float64))
    return maximum(abs.(current_values-previous_values)./denominator)
end

function response_relative_change(current::AbstractVector,previous::AbstractVector)
    length(current) == length(previous) || error(
        "Protected-response vectors changed length during coupling.",
    )
    denominator = max.(abs.(current),abs.(previous),eps(Float64))
    return maximum(abs.(current-previous)./denominator)
end

"""
    solve_closed_opensn_radiant_coupling(initial_forward, radiant_solve, opensn_solve; ...)

Iterate a two-way interface. `radiant_solve(forward_source)` returns a `Coupling_Solver_Output`
whose boundary source is the Radiant-to-OpenSn return source. `opensn_solve(return_source)` returns
an output whose boundary source is the next OpenSn-to-Radiant forward source. Linear relaxation is
applied only to the forward source. Every iterate remains nonnegative; clipping and negative-source
repair are fatal.
"""
function solve_closed_opensn_radiant_coupling(
    initial_forward::Boundary_Angular_Current_Source,
    radiant_solve::Function,
    opensn_solve::Function;
    settings::Closed_Coupling_Settings=Closed_Coupling_Settings(),
    ownership_map_hash::AbstractString,
    provenance::AbstractDict=Dict{String,String}(),
)
    isempty(ownership_map_hash) && error("Closed coupling requires an ownership-map hash.")
    forward = initial_forward
    previous_forward = nothing
    previous_return = nothing
    previous_responses = nothing
    records = Coupling_Iteration_Record[]
    final_radiant = nothing
    final_opensn = nothing
    converged = false

    for iteration in 1:settings.maximum_iterations
        radiant_output = radiant_solve(forward)
        radiant_output isa Coupling_Solver_Output || error(
            "Radiant callback must return Coupling_Solver_Output.",
        )
        opensn_output = opensn_solve(radiant_output.boundary_source)
        opensn_output isa Coupling_Solver_Output || error(
            "OpenSn callback must return Coupling_Solver_Output.",
        )
        any((radiant_output.clipping_used,opensn_output.clipping_used)) && error(
            "Closed coupling cannot use source clipping.",
        )
        any((radiant_output.negative_source_repair_used,
             opensn_output.negative_source_repair_used)) && error(
            "Closed coupling cannot repair a negative source silently.",
        )
        abs(radiant_output.particle_balance_residual) <= settings.particle_balance_tolerance ||
            error("Radiant particle balance failed during closed coupling.")
        abs(radiant_output.energy_balance_residual) <= settings.energy_balance_tolerance ||
            error("Radiant energy balance failed during closed coupling.")
        abs(opensn_output.particle_balance_residual) <= settings.particle_balance_tolerance ||
            error("OpenSn particle balance failed during closed coupling.")
        abs(opensn_output.energy_balance_residual) <= settings.energy_balance_tolerance ||
            error("OpenSn energy balance failed during closed coupling.")

        candidate_forward = opensn_output.boundary_source
        relaxed_forward = blend_boundary_sources(
            forward,candidate_forward,settings.relaxation;
            source_hash="closed-coupling-forward-$(iteration)",
        )
        responses = vcat(
            radiant_output.protected_responses,opensn_output.protected_responses,
        )
        if isnothing(previous_forward)
            forward_change = Inf
            forward_energy_change = Inf
            return_change = Inf
            return_energy_change = Inf
            protected_change = Inf
        else
            forward_change = boundary_current_relative_change(relaxed_forward,previous_forward)
            forward_energy_change = _relative_scalar_change(
                boundary_energy_current(relaxed_forward),
                boundary_energy_current(previous_forward),
            )
            return_change = boundary_current_relative_change(
                radiant_output.boundary_source,previous_return,
            )
            return_energy_change = _relative_scalar_change(
                boundary_energy_current(radiant_output.boundary_source),
                boundary_energy_current(previous_return),
            )
            protected_change = response_relative_change(responses,previous_responses)
        end
        iteration_converged = iteration >= settings.minimum_iterations &&
            forward_change <= settings.current_relative_tolerance &&
            forward_energy_change <= settings.energy_current_relative_tolerance &&
            return_change <= settings.current_relative_tolerance &&
            return_energy_change <= settings.energy_current_relative_tolerance &&
            protected_change <= settings.response_relative_tolerance
        push!(records,Coupling_Iteration_Record(
            Int64(iteration),forward_change,forward_energy_change,return_change,
            return_energy_change,protected_change,opensn_output.particle_balance_residual,
            opensn_output.energy_balance_residual,radiant_output.particle_balance_residual,
            radiant_output.energy_balance_residual,iteration_converged,
        ))
        previous_forward = relaxed_forward
        previous_return = radiant_output.boundary_source
        previous_responses = responses
        forward = relaxed_forward
        final_radiant = radiant_output
        final_opensn = opensn_output
        if iteration_converged
            converged = true
            break
        end
    end
    isnothing(final_radiant) && error("Closed coupling performed no Radiant iteration.")
    isnothing(final_opensn) && error("Closed coupling performed no OpenSn iteration.")
    provenance_string = Dict{String,String}()
    for (key,value) in provenance
        provenance_string[string(key)] = string(value)
    end
    provenance_string["schema"] = "opensn-radiant.closed-coupling/v1"
    provenance_string["relaxation"] = string(settings.relaxation)
    provenance_string["converged"] = string(converged)
    return Closed_Coupling_Result(
        forward,final_radiant.boundary_source,final_opensn,final_radiant,records,converged,
        String(ownership_map_hash),forward.normalization.source_hash,
        final_radiant.boundary_source.normalization.source_hash,provenance_string,
    )
end

function closed_coupling_receipt(result::Closed_Coupling_Result)
    last_record = last(result.iterations)
    return Dict{String,Any}(
        "schema" => "opensn-radiant.coupling/v1",
        "converged" => result.converged,
        "iteration_count" => length(result.iterations),
        "ownership_map_hash" => result.ownership_map_hash,
        "forward_source_hash" => result.forward_source_hash,
        "return_source_hash" => result.return_source_hash,
        "opensn_result_artifact_hash" => result.opensn_output.result_artifact_hash,
        "radiant_result_artifact_hash" => result.radiant_output.result_artifact_hash,
        "forward_current_closure_pass" => result.converged,
        "return_current_closure_pass" => result.converged,
        "energy_current_closure_pass" => result.converged,
        "response_convergence_pass" => result.converged,
        "clipping_used" => result.opensn_output.clipping_used ||
                           result.radiant_output.clipping_used,
        "negative_source_repair_used" =>
            result.opensn_output.negative_source_repair_used ||
            result.radiant_output.negative_source_repair_used,
        "last_forward_current_change" => last_record.forward_current_change,
        "last_forward_energy_current_change" => last_record.forward_energy_current_change,
        "last_return_current_change" => last_record.return_current_change,
        "last_return_energy_current_change" => last_record.return_energy_current_change,
        "last_protected_response_change" => last_record.protected_response_change,
        "provenance" => result.provenance,
    )
end
