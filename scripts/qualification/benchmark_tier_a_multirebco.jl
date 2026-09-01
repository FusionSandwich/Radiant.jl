#!/usr/bin/env julia

using Radiant
using SHA
using TOML

const INPUTS = [
    REBCO_Tier_A_Slab_Input(
        "YBCO",6.368105210092579,666.19,Dict("Y-89" => (1.0,1.28)),1.5,
        "78732861bcf898f84612701a64f8551e1ec41fcdb1a1692bed3a8def645b6b90",
    ),
    REBCO_Tier_A_Slab_Input(
        "GdBCO",6.942420423104715,734.54,
        Dict("Gd-155" => (0.1480,62200.0),"Gd-157" => (0.1565,239800.0)),1.0,
        "75973f9530a68559b080d7ce3e2ed692b616d26a320f41d480ca0ea28433e782",
    ),
    REBCO_Tier_A_Slab_Input(
        "EuBCO",6.900508201807441,729.25,
        Dict("Eu-151" => (0.4781,9051.0),"Eu-153" => (0.5219,364.0)),3.6,
        "70e77dbe8f2ba0e7b818fe47ee60f32c95f0438db55de67ed33dc3f41e589e56",
    ),
    REBCO_Tier_A_Slab_Input(
        "SmBCO",6.8393669104427675,727.65,Dict("Sm-149" => (0.1382,74500.0)),5.0,
        "aee8fe606e37fdfa1f9f73b52485cf1fa6b6b3a0a638997cef4fffb1b6ddb9c7",
    ),
]

median_value(values) = sort(values)[cld(length(values),2)]

function median_elapsed(operation;samples=9)
    operation()
    values = Float64[]
    for _ in 1:samples
        push!(values,@elapsed operation())
    end
    return median_value(values)
end

function minimum_equal_error_cells(input;limit=0.01,angle=1.0)
    for cells in (8,16,32,64,128,256,512,1024,2048,4096,8192)
        result = solve_rebco_tier_a_slab(
            input;angle_from_tape_plane_deg=angle,cells=cells,cache_attenuation=false,
        )
        result["peak_sublayer_capture_relative_error"] <= limit && return cells,result
    end
    error("Equal-error cell search exceeded the preregistered 8192-cell ceiling.")
end

function mapped_kernel_benchmark(input)
    target = Curved_Performance_Target(
        target_id="dpa-protected-response-tolerance-freeze-v2",maximum_cells=100_000,
    )
    material_scale = solve_rebco_tier_a_slab(
        input;angle_from_tape_plane_deg=10.0,cells=256,cache_attenuation=true,
    )["capture_fraction"]
    return benchmark_curved_pipeline(
        target,
        () -> begin
            axis_1 = collect(range(-0.02,0.02,length=5))
            axis_2 = collect(range(-0.01,0.01,length=3))
            axis_3 = collect(range(0.0,4.0*pi,length=257))
            baseline = mapped_cartesian_grid(
                axis_1,axis_2,axis_3;geometry_hash="tier-a-cartesian-$(input.material_tag)",
            )
            improved = mapped_helical_grid(
                axis_1,axis_2,axis_3;radius_cm=2.0,pitch_per_turn_cm=0.5,
                coordinate_system=:repeated_turns,
                geometry_hash="tier-a-curved-$(input.material_tag)",
            )
            (baseline=baseline,improved=improved)
        end;
        baseline_sweep=grid -> mapped_streaming_divergence(
            grid,[0.0,0.0,1.0],fill(material_scale,size(grid.cell_volumes_cm3));
            boundary_inflow=(axis,side,i,j,k) -> material_scale,
        ),
        improved_sweep=grid -> mapped_streaming_divergence(
            grid,[0.0,0.0,1.0],fill(material_scale,size(grid.cell_volumes_cm3));
            boundary_inflow=(axis,side,i,j,k) -> material_scale,
        ),
        scoring=(result,grid) -> sum(abs,result.angular_divergence),
        baseline_remap=(result,grid) -> reverse(reverse(result.angular_divergence,dims=3),dims=3),
        cached_atlas_remap=(result,grid) -> reverse(result.angular_divergence,dims=3),
        samples=9,
    )
end

function family_report(input)
    equal_cells,equal_result = minimum_equal_error_cells(input)
    equal_cartesian_s = median_elapsed(() -> solve_rebco_tier_a_slab(
        input;angle_from_tape_plane_deg=1.0,cells=equal_cells,cache_attenuation=false,
    ))
    equal_cached_s = median_elapsed(() -> solve_rebco_tier_a_slab(
        input;angle_from_tape_plane_deg=1.0,cells=equal_cells,cache_attenuation=true,
    ))
    mapped = mapped_kernel_benchmark(input)
    angular = Dict{String,Any}()
    for angle in (90.0,45.0,10.0,5.0,1.0)
        result = solve_rebco_tier_a_slab(
            input;angle_from_tape_plane_deg=angle,cells=512,cache_attenuation=true,
        )
        angular[string(angle)] = Dict(
            "capture_fraction" => result["capture_fraction"],
            "integrated_capture_relative_error" => result["integrated_capture_relative_error"],
            "peak_sublayer_capture_relative_error" =>
                result["peak_sublayer_capture_relative_error"],
        )
    end
    return Dict{String,Any}(
        "material_packet_sha256" => input.material_packet_sha256,
        "film_thickness_um" => input.film_thickness_um,
        "macroscopic_capture_cm_inv" => rebco_tier_a_macroscopic_capture(input),
        "source_energy_eV" => input.source_energy_eV,
        "mapped_curved_kernel" => mapped,
        "equal_error" => Dict(
            "protected_response" => "peak sublayer rare-earth capture density",
            "relative_error_limit" => 0.01,"cells" => equal_cells,
            "protected_response_relative_error" =>
                equal_result["peak_sublayer_capture_relative_error"],
            "cartesian_s" => equal_cartesian_s,"cached_local_s" => equal_cached_s,
            "cached_to_cartesian_ratio" => equal_cached_s/max(equal_cartesian_s,eps(Float64)),
            "iteration_count" => 1,
        ),
        "angular_capture" => angular,
        "reference_physical_candidate" => true,"lot_specific" => false,
        "production_qualified" => false,
    )
end

output = length(ARGS) == 1 ? abspath(ARGS[1]) :
    joinpath(pwd(),"qualification","TIER_A_MULTIREBCO_PERFORMANCE_RAW.toml")
families = Dict(input.material_tag => family_report(input) for input in INPUTS)
report = Dict{String,Any}(
    "schema" => "radiant.hts.tier_a_multirebco_performance_raw/v1",
    "julia_version" => string(VERSION),"radiant_project" => abspath(pwd()),
    "evidence_tier" => "A_REFERENCE_PHYSICAL_CANDIDATE",
    "protected_tolerance_freeze_id" => "RADIANTHTS_PROTECTED_RESPONSE_TOLERANCE_FREEZE_V2",
    "preserved_manufactured_curved_cartesian_ratio" => 1.1133,
    "preregistered_equal_unknown_limit" => 1.25,
    "preregistered_equal_error_limit" => 1.35,
    "candidate_thickness_axis_um" => [1.0,0.5,0.25,0.1],
    "source_faces" => ["front","rear","edge_1","edge_2","end_1","end_2"],
    "face_coverage_status" => "LOCAL_SLAB_SYMMETRY_ONLY_NOT_FINITE_TAPE_TRANSPORT",
    "families" => families,
    "physical_transport_qualification" => false,
)
mkpath(dirname(output))
open(output,"w") do io
    TOML.print(value -> value === nothing ? "NOT_MEASURED" : value,io,report;sorted=true)
end
println("report=$(output)")
println("sha256=$(bytes2hex(sha256(read(output))))")
