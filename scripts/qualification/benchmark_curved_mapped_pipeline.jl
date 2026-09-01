#!/usr/bin/env julia

using Radiant

# Frozen before any timed work. This is a bounded software-performance target, not a physical gate.
target = Curved_Performance_Target(
    target_id="dpa-b653f509-curved-performance-target",
    maximum_cells=100_000,
)

result = benchmark_curved_pipeline(
    target,
    () -> begin
        axis_1 = collect(range(-0.02,0.02,length=5))
        axis_2 = collect(range(-0.01,0.01,length=3))
        axis_3 = collect(range(0.0,4.0*pi,length=257))
        baseline = mapped_cartesian_grid(
            axis_1,axis_2,axis_3;geometry_hash="matched-cartesian-performance-fixture",
        )
        improved = mapped_helical_grid(
            axis_1,axis_2,axis_3;radius_cm=2.0,pitch_per_turn_cm=0.5,
            coordinate_system=:repeated_turns,
            geometry_hash="analytic-curved-performance-fixture",
        )
        (baseline=baseline,improved=improved)
    end;
    baseline_sweep=grid -> mapped_streaming_divergence(
        grid,[0.0,0.0,1.0],ones(size(grid.cell_volumes_cm3));
        boundary_inflow=(axis,side,i,j,k) -> 1.0,
    ),
    improved_sweep=grid -> mapped_streaming_divergence(
        grid,[0.0,0.0,1.0],ones(size(grid.cell_volumes_cm3));
        boundary_inflow=(axis,side,i,j,k) -> 1.0,
    ),
    scoring=(result,grid) -> sum(abs,result.angular_divergence),
    baseline_remap=(result,grid) ->
        reverse(reverse(result.angular_divergence,dims=3),dims=3),
    cached_atlas_remap=(result,grid) -> reverse(result.angular_divergence,dims=3),
)

for key in sort(collect(keys(result)))
    println("$(key)=$(result[key])")
end
result["measurement_status"]["equal_error_solve"] == "NOT_MEASURED_DIAGNOSTIC_ONLY" ||
    error("Equal-error solve must remain unmeasured until a protected-response solve is added.")
result["passed"] || error(
    "Curved mapped benchmark exceeded a measured kernel, memory, or remap target.",
)
