using Radiant
using Test

@testset "Hash-bound HDF5 source interchange" begin
    photon = Photon()
    normalization = Source_Normalization(
        basis=:per_history,
        source_rate_per_s=2.0,
        symmetry_factor=4.0,
        time_interval_s=(0.0,10.0),
        time_class=:prompt,
        source_hash="hdf5-roundtrip-source",
        provenance=Dict("global_model_hash" => "global-fixture"),
    )

    @testset "Boundary phase-space round trip" begin
        patch_ids = [1,2,3]
        centroids = [
            0.0 0.0 0.0
            0.0 1.0 0.0
            0.0 2.0 0.0
        ]
        areas = [1.0,2.0,3.0]
        normals = repeat([1.0 0.0 0.0],3,1)
        tangent_1 = repeat([0.0 1.0 0.0],3,1)
        tangent_2 = repeat([0.0 0.0 1.0],3,1)
        edges = [1.0,10.0,100.0]
        directions = [
            -1.0 0.0 0.0
            -0.5 sqrt(3.0)/2.0 0.0
        ]
        weights = [0.25,0.75]
        angular_flux = zeros(Float64,3,2,2)
        for patch in 1:3, group in 1:2
            angular_flux[patch,group,:] .= [patch+group,2.0*patch+group]
        end
        variance = 0.01 .* angular_flux.^2
        source = Boundary_Angular_Current_Source(
            photon,patch_ids,centroids,areas,normals,tangent_1,tangent_2,
            edges,directions,weights,angular_flux,normalization;
            variance=variance,
            provenance=Dict("producer" => "roundtrip-test"),
        )

        mktempdir() do directory
            path = joinpath(directory,"boundary.h5")
            artifact_hash = write_boundary_angular_current_hdf5(path,source)
            @test length(artifact_hash) == 64
            restored = read_boundary_angular_current_hdf5(
                path,photon;expected_file_sha256=artifact_hash,
            )
            @test restored.patch_ids == source.patch_ids
            @test restored.centroids_cm == source.centroids_cm
            @test restored.normals == source.normals
            @test restored.tangent_1 == source.tangent_1
            @test restored.tangent_2 == source.tangent_2
            @test restored.energy_edges_eV == source.energy_edges_eV
            @test restored.directions == source.directions
            @test restored.quadrature_weights == source.quadrature_weights
            @test restored.angular_flux == source.angular_flux
            @test restored.variance == source.variance
            @test get_incoming_current(restored) ≈ get_incoming_current(source)
            @test restored.provenance["file_sha256"] == artifact_hash
            @test_throws ErrorException read_boundary_angular_current_hdf5(
                path,photon;expected_file_sha256=repeat("0",64),
            )
        end
    end

    @testset "Anisotropic volume source round trip" begin
        values = zeros(Float64,3,2,2)
        for voxel in 1:3, group in 1:2
            values[voxel,group,:] .= [voxel+group,0.5*(voxel+group)]
        end
        directions = [
            -1.0 0.0 0.0
            1.0 0.0 0.0
        ]
        weights = [1.0,1.0]
        source = Anisotropic_Volume_Source(
            photon,[1,2,3],[1.0,2.0,3.0],[1.0,10.0,100.0],:ordinates,
            values,normalization;
            directions=directions,
            quadrature_weights=weights,
            variance=0.02 .* values.^2,
            parent_reaction="synthetic-neutron-to-photon",
            provenance=Dict("producer" => "roundtrip-test"),
        )

        mktempdir() do directory
            path = joinpath(directory,"volume.h5")
            artifact_hash = write_anisotropic_volume_source_hdf5(path,source)
            restored = read_anisotropic_volume_source_hdf5(
                path,photon;expected_file_sha256=artifact_hash,
            )
            @test restored.voxel_ids == source.voxel_ids
            @test restored.voxel_volumes_cm3 == source.voxel_volumes_cm3
            @test restored.energy_edges_eV == source.energy_edges_eV
            @test restored.angular_representation == :ordinates
            @test restored.directions == source.directions
            @test restored.quadrature_weights == source.quadrature_weights
            @test restored.values == source.values
            @test restored.variance == source.variance
            @test restored.parent_reaction == source.parent_reaction
            @test get_volume_source_rate(restored) ≈ get_volume_source_rate(source)
            @test restored.provenance["file_sha256"] == artifact_hash
        end
    end
end
