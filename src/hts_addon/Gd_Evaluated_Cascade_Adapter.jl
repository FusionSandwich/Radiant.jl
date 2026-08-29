"""Normalized branch-level Gd cascade data and its hash-bound lineage."""
struct Evaluated_Gd_Cascade_Dataset
    schema::String
    cascades::Dict{String,Gd_Prompt_Capture_Cascade}
    branch_probabilities::Dict{String,Dict{String,Float64}}
    input_file_sha256::String
    evaluation_id::String
    evaluation_version::String
    provenance::Dict{String,String}
end

function _cascade_species(value::AbstractString)
    normalized = lowercase(replace(strip(value),'-' => '_',' ' => '_'))
    mapping = Dict(
        "photon" => :photon,
        "gamma" => :photon,
        "xray" => :xray,
        "x_ray" => :xray,
        "conversion_electron" => :conversion_electron,
        "internal_conversion_electron" => :conversion_electron,
        "auger_electron" => :auger_electron,
    )
    haskey(mapping,normalized) || error("Unsupported evaluated cascade species $(value).")
    return mapping[normalized]
end

"""
    read_evaluated_gd_cascade_csv(path; ...)

Read a normalized branch-resolved Gd-155/Gd-157 capture cascade. Required columns are:

```text
target_nuclide,residual_nuclide,branch_id,branch_probability,q_value_eV,
species,energy_eV,yield_per_branch,delay_s,correlation_id,transition_id,
recoil_energy_eV
```

Each branch probability set must sum to one and each branch must close its Q value after all
photons, x-rays, conversion electrons, Auger electrons, and recoil are included. Branch-weighted
emission yields are converted to the mean-per-capture source representation used by Radiant.
Incomplete gamma-only data therefore fail rather than sending missing atomic/recoil energy to heat.
"""
function read_evaluated_gd_cascade_csv(
    path::AbstractString;
    evaluation_id::AbstractString,
    evaluation_version::AbstractString,
    expected_file_sha256::Union{Nothing,AbstractString}=nothing,
    qualification_status::Symbol=:candidate,
    probability_tolerance::Real=1.0e-8,
    energy_rtol::Real=5.0e-3,
    energy_atol_eV::Real=1.0e3,
    provenance::AbstractDict=Dict{String,String}(),
)
    isempty(evaluation_id) && error("Gd evaluation identifier cannot be empty.")
    isempty(evaluation_version) && error("Gd evaluation version cannot be empty.")
    artifact_hash = _verify_file_hash(path,expected_file_sha256)
    header,rows = _simple_csv_rows(path)
    required = [
        "target_nuclide","residual_nuclide","branch_id","branch_probability",
        "q_value_eV","species","energy_eV","yield_per_branch","delay_s",
        "correlation_id","transition_id","recoil_energy_eV",
    ]
    all(key in header for key in required) || error(
        "Evaluated Gd cascade CSV is missing required columns.",
    )
    targets = sort(unique(row["target_nuclide"] for row in rows))
    all(target in ("Gd-155","Gd-157") for target in targets) || error(
        "Evaluated Gd cascade adapter accepts only isotope-resolved Gd-155 or Gd-157 targets.",
    )

    cascades = Dict{String,Gd_Prompt_Capture_Cascade}()
    all_probabilities = Dict{String,Dict{String,Float64}}()
    provenance_string = Dict{String,String}(
        "source_path" => abspath(path),
        "source_sha256" => artifact_hash,
        "evaluation_id" => String(evaluation_id),
        "evaluation_version" => String(evaluation_version),
        "branch_energy_closure" => "required",
        "atomic_relaxation_included" => "required",
    )
    for (key,value) in provenance
        provenance_string[string(key)] = string(value)
    end

    for target in targets
        target_rows = [row for row in rows if row["target_nuclide"] == target]
        residuals = unique(row["residual_nuclide"] for row in target_rows)
        length(residuals) == 1 || error("Target $(target) maps to multiple residual nuclides.")
        q_values = unique(_csv_float(row,"q_value_eV") for row in target_rows)
        length(q_values) == 1 || error("Target $(target) has inconsistent Q values.")
        q_value = first(q_values)
        branches = sort(unique(row["branch_id"] for row in target_rows))
        probabilities = Dict{String,Float64}()
        weighted_emissions = Capture_Emission_Line[]
        mean_recoil = 0.0

        for branch_id in branches
            branch_rows = [row for row in target_rows if row["branch_id"] == branch_id]
            probability_values = unique(_csv_float(row,"branch_probability")
                                        for row in branch_rows)
            length(probability_values) == 1 || error(
                "Branch $(branch_id) has inconsistent probabilities.",
            )
            probability = first(probability_values)
            0.0 <= probability <= 1.0 || error("Gd branch probability lies outside [0,1].")
            probabilities[branch_id] = probability
            recoil_values = unique(_csv_float(row,"recoil_energy_eV") for row in branch_rows)
            length(recoil_values) == 1 || error(
                "Branch $(branch_id) has inconsistent recoil energy.",
            )
            recoil = first(recoil_values)
            recoil >= 0.0 || error("Gd branch recoil energy cannot be negative.")
            branch_emitted_energy = 0.0
            for row in branch_rows
                species = _cascade_species(row["species"])
                energy = _csv_float(row,"energy_eV")
                yield_value = _csv_float(row,"yield_per_branch")
                delay = _csv_float(row,"delay_s")
                yield_value >= 0.0 || error("Gd emission yield cannot be negative.")
                delay >= 0.0 || error("Gd prompt-emission delay cannot be negative.")
                branch_emitted_energy += energy*yield_value
                push!(weighted_emissions,Capture_Emission_Line(
                    species,energy,probability*yield_value;
                    delay_s=delay,
                    correlation_id=string(branch_id,"/",row["correlation_id"]),
                    transition_id=row["transition_id"],
                    metadata=Dict(
                        "branch_id" => branch_id,
                        "branch_probability" => string(probability),
                        "yield_per_branch" => string(yield_value),
                        "evaluation_id" => String(evaluation_id),
                    ),
                ))
            end
            isapprox(branch_emitted_energy+recoil,q_value;
                     rtol=energy_rtol,atol=energy_atol_eV) || error(
                "Gd branch $(branch_id) fails Q-value closure: emitted=$(branch_emitted_energy), recoil=$(recoil), Q=$(q_value).",
            )
            mean_recoil += probability*recoil
        end
        isapprox(sum(values(probabilities)),1.0;
                 rtol=probability_tolerance,atol=probability_tolerance) || error(
            "Gd branch probabilities for $(target) do not sum to one.",
        )
        all_probabilities[target] = probabilities
        cascade = Gd_Prompt_Capture_Cascade(
            target,first(residuals),q_value,weighted_emissions;
            mean_recoil_energy_eV=mean_recoil,
            evaluation_id=evaluation_id,
            evaluation_hash=artifact_hash,
            qualification_status=qualification_status,
            metadata=merge(copy(provenance_string),Dict(
                "target_nuclide" => target,
                "branch_count" => string(length(branches)),
                "mean_per_capture" => "true",
            )),
        )
        emitted_mean = sum(line.energy_eV*line.yield_per_capture
                           for line in cascade.emissions)
        isapprox(emitted_mean+mean_recoil,q_value;
                 rtol=energy_rtol,atol=energy_atol_eV) || error(
            "Branch-weighted cascade for $(target) fails mean Q-value closure.",
        )
        cascades[target] = cascade
    end
    return Evaluated_Gd_Cascade_Dataset(
        "radiant.hts.evaluated_gd_cascade/v1",cascades,all_probabilities,
        artifact_hash,String(evaluation_id),String(evaluation_version),provenance_string,
    )
end

function evaluated_gd_dataset_is_production_ready(dataset::Evaluated_Gd_Cascade_Dataset)
    required = Set(["Gd-155","Gd-157"])
    return Set(keys(dataset.cascades)) == required &&
           all(cascade.qualification_status == :qualified
               for cascade in values(dataset.cascades)) &&
           length(dataset.input_file_sha256) == 64 &&
           get(dataset.provenance,"source_classification","") in
               ("evaluated-nuclear-data","measured-cascade-repository")
end

function gd_dataset_receipt(dataset::Evaluated_Gd_Cascade_Dataset)
    isotopes = Dict{String,Any}()
    for (target,cascade) in dataset.cascades
        isotopes[target] = Dict{String,Any}(
            "residual_nuclide" => cascade.residual_nuclide,
            "q_value_eV" => cascade.q_value_eV,
            "mean_recoil_energy_eV" => cascade.mean_recoil_energy_eV,
            "emission_count" => length(cascade.emissions),
            "branch_probabilities" => dataset.branch_probabilities[target],
            "qualification_status" => string(cascade.qualification_status),
        )
    end
    return Dict{String,Any}(
        "schema" => dataset.schema,
        "evaluation_id" => dataset.evaluation_id,
        "evaluation_version" => dataset.evaluation_version,
        "input_file_sha256" => dataset.input_file_sha256,
        "isotopes" => isotopes,
        "production_ready" => evaluated_gd_dataset_is_production_ready(dataset),
        "provenance" => copy(dataset.provenance),
    )
end

export Evaluated_Gd_Cascade_Dataset,read_evaluated_gd_cascade_csv
export evaluated_gd_dataset_is_production_ready,gd_dataset_receipt
