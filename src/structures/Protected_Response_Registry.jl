"""One response whose definition and convergence criteria are frozen before a calculation."""
struct Protected_Response
    id::String
    domain::String
    particle::String
    time_state::String
    production_owner::String
    units::String
    normalization::String
    angular_representation::String
    reference_solver::String
    relative_tolerance::Float64
    absolute_tolerance::Float64
    uncertainty_components::Vector{String}

    function Protected_Response(;
        id::AbstractString,
        domain::AbstractString,
        particle::AbstractString,
        time_state::AbstractString="prompt",
        production_owner::AbstractString,
        units::AbstractString,
        normalization::AbstractString,
        angular_representation::AbstractString="scalar",
        reference_solver::AbstractString,
        relative_tolerance::Real,
        absolute_tolerance::Real=0.0,
        uncertainty_components::AbstractVector{<:AbstractString}=String[],
    )
        isempty(id) && error("Protected-response identifier cannot be empty.")
        isempty(domain) && error("Protected-response domain cannot be empty.")
        isempty(particle) && error("Protected-response particle cannot be empty.")
        isempty(production_owner) && error("Protected response requires a production owner.")
        isempty(units) && error("Protected-response units cannot be empty.")
        isempty(normalization) && error("Protected-response normalization cannot be empty.")
        isempty(reference_solver) && error("Protected response requires a reference solver.")
        relative = Float64(relative_tolerance)
        absolute = Float64(absolute_tolerance)
        isfinite(relative) && relative ≥ 0.0 || error(
            "Protected-response relative tolerance must be finite and nonnegative.",
        )
        isfinite(absolute) && absolute ≥ 0.0 || error(
            "Protected-response absolute tolerance must be finite and nonnegative.",
        )
        return new(
            String(id),String(domain),String(particle),String(time_state),
            String(production_owner),String(units),String(normalization),
            String(angular_representation),String(reference_solver),relative,absolute,
            String[string(value) for value in uncertainty_components],
        )
    end
end

struct Protected_Response_Registry
    schema::String
    responses::Vector{Protected_Response}

    function Protected_Response_Registry(
        responses::AbstractVector{Protected_Response};
        schema::AbstractString="radiant.hts.protected_responses/v1",
    )
        isempty(responses) && error("Protected-response registry cannot be empty.")
        response_vector = Protected_Response[responses...]
        ids = getfield.(response_vector,:id)
        length(unique(ids)) == length(ids) || error(
            "Protected-response identifiers must be unique.",
        )
        return new(String(schema),response_vector)
    end
end

function get_protected_response(this::Protected_Response_Registry,id::AbstractString)
    index = findfirst(response -> response.id == id,this.responses)
    isnothing(index) && error("No protected response matches identifier $(id).")
    return this.responses[index]
end

function response_is_converged(
    response::Protected_Response,
    reference::Real,
    candidate::Real,
)
    isfinite(reference) && isfinite(candidate) || return false
    return isapprox(
        Float64(candidate),Float64(reference);
        rtol=response.relative_tolerance,
        atol=response.absolute_tolerance,
    )
end

"""Core HTS, fibre, and detector responses used by the local deterministic qualification."""
function default_hts_protected_responses()
    common_uncertainty = [
        "source-bank sampling",
        "source projection",
        "energy grouping",
        "angular quadrature",
        "spatial mesh",
        "cross-section data",
        "normalization",
    ]
    responses = Protected_Response[
        Protected_Response(
            id="ybco_em_deposition",
            domain="YBCO layer",
            particle="photon+electron+positron",
            production_owner="Radiant",
            units="eV/source history",
            normalization="explicit Source_Normalization",
            reference_solver="Geant4 matched-physics",
            relative_tolerance=0.05,
            uncertainty_components=common_uncertainty,
        ),
        Protected_Response(
            id="ybco_peak_prompt_heat",
            domain="YBCO sublayer voxel",
            particle="all owned prompt particles",
            production_owner="ownership-map aggregate",
            units="W/cm^3",
            normalization="full-device source rate",
            reference_solver="local continuous-energy OpenMC/Geant4",
            relative_tolerance=0.10,
            uncertainty_components=vcat(common_uncertainty,["energy partition model"]),
        ),
        Protected_Response(
            id="ag_to_ybco_electron_current",
            domain="Ag/YBCO interface",
            particle="electron",
            production_owner="Radiant",
            units="electrons/source history",
            normalization="interface-current integral",
            angular_representation="outgoing half-range current",
            reference_solver="Geant4 matched-physics",
            relative_tolerance=0.05,
            uncertainty_components=common_uncertainty,
        ),
        Protected_Response(
            id="fibre_core_ionizing_dose",
            domain="fibre core",
            particle="photon+electron+positron",
            production_owner="Radiant",
            units="Gy",
            normalization="mass and physical source history",
            reference_solver="Geant4 matched-physics",
            relative_tolerance=0.05,
            uncertainty_components=common_uncertainty,
        ),
        Protected_Response(
            id="fibre_capture_and_displacement_source",
            domain="fibre assembly",
            particle="neutron",
            production_owner="OpenSn/OpenMC",
            units="reactions/source history",
            normalization="volume scalar flux and nuclide density",
            reference_solver="continuous-energy OpenMC",
            relative_tolerance=0.05,
            uncertainty_components=vcat(common_uncertainty,["thermal scattering","self shielding"]),
        ),
        Protected_Response(
            id="hts_detector_active_energy_impulse",
            domain="detector active HTS region",
            particle="all converter and EM products",
            production_owner="Radiant+external reaction owner without overlap",
            units="eV/event",
            normalization="weighted event kernel",
            reference_solver="Geant4 event-wise reference",
            relative_tolerance=0.05,
            uncertainty_components=vcat(common_uncertainty,["event-kernel sampling"]),
        ),
        Protected_Response(
            id="hts_detector_pulse_height",
            domain="detector electrothermal model",
            particle="response state",
            production_owner="electrothermal response solver",
            units="V or A",
            normalization="per accepted detector event",
            reference_solver="calibrated transient experiment/model",
            relative_tolerance=0.10,
            uncertainty_components=[
                "energy-deposition impulse",
                "heat capacity",
                "interface conductance",
                "transition curve",
                "bias state",
                "readout model",
            ],
        ),
    ]
    return Protected_Response_Registry(responses)
end
