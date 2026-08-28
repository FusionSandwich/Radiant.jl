using Test
using LinearAlgebra

@testset "Hash-bound HDF5 source interchange" begin
    temporary = mktempdir()
    normalization = Source_Normalization(
        basis=:per_history,
        source_rate_per_s=3.0,
        symmetry_factor=2.0,
        time_interval_s=(0.0,10.0),
        time_class=:prompt,
        source_hash="producer-boundary-hash",
        provenance=Dict("producer" => "analytic-openmc-fixture"),
    )

    particle = Photon()
    patch_ids = [7,9]
    centroids = [0.0 0.0 0.0; 1.0 0.0 0.0]
    areas = [1.0,2.0]
    normals = [-1.0 0.0 0.0; -1.0 0.0 0.0]
    tangent_1 = [0.0 1.0 0.0; 0.0 1.0 0.0]
    tangent_2 = [0.0 0.0 -1.0; 0.0 0.0 -1.0]
    energy_edges = [1.0e3,1.0e4,1.0e6]
    directions = [1.0 0.0 0.0; 0.5 sqrt(3.0)/2.0 0.0]
    weights = [0.25,0.75]
    angular_flux = reshape(collect(1.0:8.0),2,2,2)
    variance = 0.01 .* angular_flux

    boundary = Boundary_Angular_Current_Source(
        particle,patch_ids,centroids,areas,normals,tangent_1,tangent_2,
        energy_edges,directions,weights,angular_flux,normalization;
        variance=variance,
        provenance=Dict("surface_semantics" => "incoming-current"),
    )
    boundary_path = joinpath(temporary,"boundary_source.h5")
    boundary_digest = write_boundary_source_hdf5(boundary_path,boundary)
    @test boundary_digest == source_artifact_sha256(boundary_path)
    replayed_boundary = read_boundary_source_hdf5(
        boundary_path;expected_sha256=boundary_digest,
    )
    @test replayed_boundary.patch_ids == boundary.patch_ids
    @test replayed_boundary.centroids_cm == boundary.centroids_cm
    @test replayed_boundary.normals == boundary.normals
    @test replayed_boundary.energy_edges_eV == boundary.energy_edges_eV
    @test replayed_boundary.directions == boundary.directions
    @test replayed_boundary.quadrature_weights == boundary.quadrature_weights
    @test replayed_boundary.angular_flux == boundary.angular_flux
    @test replayed_boundary.variance == boundary.variance
    @test get_incoming_current(replayed_boundary) ≈ get_incoming_current(boundary)
    @test replayed_boundary.provenance["artifact_sha256"] == boundary_digest
    @test_throws ErrorException read_boundary_source_hdf5(
        boundary_path;expected_sha256=repeat("0",64),
    )
    @test read_radiant_source_hdf5(boundary_path) isa Boundary_Angular_Current_Source

    volume_normalization = Source_Normalization(
        basis=:per_second,
        source_rate_per_s=1.0,
        symmetry_factor=1.0,
        time_interval_s=(100.0,200.0),
        time_class=:delayed,
        source_hash="activation-source-hash",
    )
    volume_values = reshape(collect(1.0:12.0),2,2,3)
    volume_variance = 0.02 .* volume_values
    volume = Anisotropic_Volume_Source(
        particle,
        [1,4],
        [0.5,1.5],
        energy_edges,
        :ordinates,
        volume_values,
        volume_normalization;
        directions=[
            1.0 0.0 0.0
            0.0 1.0 0.0
            0.0 0.0 1.0
        ],
        quadrature_weights=[0.2,0.3,0.5],
        variance=volume_variance,
        parent_reaction="decay-photon",
        provenance=Dict("cooling_time_s" => "100.0"),
    )
    volume_path = joinpath(temporary,"volume_source.h5")
    volume_digest = write_volume_source_hdf5(volume_path,volume)
    replayed_volume = read_volume_source_hdf5(
        volume_path;expected_sha256=volume_digest,
    )
    @test replayed_volume.voxel_ids == volume.voxel_ids
    @test replayed_volume.voxel_volumes_cm3 == volume.voxel_volumes_cm3
    @test replayed_volume.angular_representation == :ordinates
    @test replayed_volume.values == volume.values
    @test replayed_volume.variance == volume.variance
    @test replayed_volume.directions == volume.directions
    @test replayed_volume.quadrature_weights == volume.quadrature_weights
    @test replayed_volume.parent_reaction == "decay-photon"
    @test replayed_volume.normalization.time_class == :delayed
    @test replayed_volume.provenance["artifact_sha256"] == volume_digest
    @test get_volume_source_rate(replayed_volume) ≈ get_volume_source_rate(volume)
    @test read_radiant_source_hdf5(volume_path) isa Anisotropic_Volume_Source
end
