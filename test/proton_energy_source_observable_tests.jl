using Test

function _energy_observable_coefficients(voxel,group,angular)
    a = 1.0 + 0.13voxel + 0.07group + 0.03angular
    b = (-1.0)^(voxel+group+angular) * (0.18 + 0.02angular)
    c = 0.04 + 0.01voxel
    return a,b,c
end

function _energy_observable_fixture(representation; energy_modes=2, physical_basis=:per_source_particle,
                                   source_rate=5.0, symmetry_factor=2.0,
                                   provenance=Dict{String,String}(), angular_count=1,
                                   weights=Float64[])
    edges = [0.7e6,3.2e6,8.9e6]
    volumes = [0.4,1.7]
    nangular = representation == :isotropic ? 1 : angular_count
    values = zeros(Float64,2,2,nangular)
    moments = zeros(Float64,2,2,nangular,energy_modes)
    for voxel in 1:2, group in 1:2, angular in 1:nangular
        a,b,c = _energy_observable_coefficients(voxel,group,angular)
        width = edges[group+1]-edges[group]
        q0 = width*(a+c/3)
        values[voxel,group,angular] = q0
        moments[voxel,group,angular,1] = q0
        if energy_modes >= 2
            moments[voxel,group,angular,2] = width*b/sqrt(3.0)
        end
    end
    normalization = Source_Normalization(basis=physical_basis,source_rate_per_s=source_rate,
        symmetry_factor=symmetry_factor,source_hash="energy-observable-test")
    directions = representation == :ordinates ? [1.0 0.0 0.0; -1.0 0.0 0.0][1:nangular,:] : zeros(Float64,0,3)
    qweights = representation == :ordinates ? weights : Float64[]
    base = Anisotropic_Volume_Source(Photon(),[2,7],volumes,edges,representation,values,
        normalization;directions=directions,quadrature_weights=qweights,provenance=provenance,
        source_tolerance=0.0)
    return Energy_Moment_Volume_Source(base,moments),edges,volumes
end

function _energy_observable_simpson(a,b,c,lower,upper)
    width = upper-lower
    midpoint = lower/2+upper/2
    q(u) = a+b*u+c*u^2
    # Three-point Simpson rule integrates E*q(E) exactly for this quadratic q.
    width/6 * (lower*q(-1.0) + 4midpoint*q(0.0) + upper*q(1.0))
end

function _energy_observable_expected(representation,edges,volumes;nangular=1,weights=Float64[],zeroth_index=1)
    result = zeros(Float64,length(edges)-1)
    for group in eachindex(result), voxel in eachindex(volumes)
        for angular in 1:nangular
            representation == :moments && angular != zeroth_index && continue
            a,b,c = _energy_observable_coefficients(voxel,group,angular)
            angular_weight = representation == :ordinates ? weights[angular] : 1.0
            result[group] += volumes[voxel]*angular_weight*
                _energy_observable_simpson(a,b,c,edges[group],edges[group+1])
        end
    end
    return result .* 1.0e-6
end

@testset "Declared volume-source first energy moment" begin
    for representation in (:isotropic,:ordinates,:moments)
        if representation == :ordinates
            weights = [0.25,1.75]
            provenance = Dict{String,String}()
            angular_count = 2
        elseif representation == :moments
            weights = Float64[]
            angular_count = 2
            provenance = Dict("angular_basis"=>"test-angle-basis",
                "zeroth_moment_index"=>"2","zeroth_moment_is_angle_integrated"=>"true")
        else
            weights = Float64[]
            angular_count = 1
            provenance = Dict{String,String}()
        end
        wrapper,edges,volumes = _energy_observable_fixture(representation;
            provenance=provenance,angular_count=angular_count,weights=weights)
        expected = _energy_observable_expected(representation,edges,volumes;
            nangular=angular_count,weights=weights,zeroth_index=representation == :moments ? 2 : 1)
        observed = Radiant.get_volume_source_energy_MeV(wrapper)
        @test observed ≈ expected rtol=2e-14 atol=2e-14
        @test observed[1] != observed[2]
        if representation == :ordinates
            unweighted = _energy_observable_expected(representation,edges,volumes;
                nangular=angular_count,weights=[1.0,1.0])
            @test observed != unweighted
        end
        @test Radiant.get_volume_source_energy_MeV(wrapper;physical=true) ≈ 10 .* expected
    end

    # A one-mode source has no Q1; its declared within-group energy is the midpoint.
    wrapper,edges,volumes = _energy_observable_fixture(:isotropic;energy_modes=1)
    midpoint_expected = [sum(volumes[v] *
        (edges[g]/2+edges[g+1]/2) * wrapper.energy_moments[v,g,1,1]
        for v in eachindex(volumes)) / 1.0e6 for g in 1:2]
    @test Radiant.get_volume_source_energy_MeV(wrapper) ≈ midpoint_expected

    # Already-per-second values apply symmetry, but never multiply by source rate again.
    per_second,edges,volumes = _energy_observable_fixture(:isotropic;physical_basis=:per_second,
        source_rate=5.0,symmetry_factor=3.0)
    expected = _energy_observable_expected(:isotropic,edges,volumes)
    @test Radiant.get_volume_source_energy_MeV(per_second;physical=true) ≈ 3 .* expected
end

@testset "Energy source observable validates angular metadata and mutable contents" begin
    edges = [0.7e6,3.2e6,8.9e6]
    volumes = [0.4,1.7]
    moments_source = Dict("angular_basis"=>"test-angle-basis",
        "zeroth_moment_index"=>"2","zeroth_moment_is_angle_integrated"=>"true")
    wrapper,_,_ = _energy_observable_fixture(:moments;provenance=moments_source,angular_count=2)
    @test isfinite(sum(Radiant.get_volume_source_energy_MeV(wrapper)))

    missing_angle = Dict("angular_basis"=>"test-angle-basis","zeroth_moment_index"=>"2")
    bad,_,_ = _energy_observable_fixture(:moments;provenance=missing_angle,angular_count=2)
    @test_throws ErrorException Radiant.get_volume_source_energy_MeV(bad)
    bad_index = Dict("angular_basis"=>"test-angle-basis","zeroth_moment_index"=>"3",
        "zeroth_moment_is_angle_integrated"=>"true")
    bad,_,_ = _energy_observable_fixture(:moments;provenance=bad_index,angular_count=2)
    @test_throws ErrorException Radiant.get_volume_source_energy_MeV(bad)
    malformed_index = Dict("angular_basis"=>"test-angle-basis","zeroth_moment_index"=>"bad",
        "zeroth_moment_is_angle_integrated"=>"true")
    bad,_,_ = _energy_observable_fixture(:moments;provenance=malformed_index,angular_count=2)
    @test_throws ErrorException Radiant.get_volume_source_energy_MeV(bad)

    bad,_,_ = _energy_observable_fixture(:isotropic)
    bad.source.energy_edges_eV[2] = bad.source.energy_edges_eV[1]
    @test_throws ErrorException Radiant.get_volume_source_energy_MeV(bad)
    bad,_,_ = _energy_observable_fixture(:isotropic)
    resize!(bad.source.energy_edges_eV,2)
    @test_throws ErrorException Radiant.get_volume_source_energy_MeV(bad)
    bad,_,_ = _energy_observable_fixture(:isotropic)
    bad.source.values[1,1,1] += 1.0
    @test_throws ErrorException Radiant.get_volume_source_energy_MeV(bad)
    bad,_,_ = _energy_observable_fixture(:isotropic)
    bad.energy_moments[1,1,1,2] = Inf
    @test_throws ErrorException Radiant.get_volume_source_energy_MeV(bad)
    bad,_,_ = _energy_observable_fixture(:isotropic)
    bad.source.voxel_volumes_cm3[1] = NaN
    @test_throws ErrorException Radiant.get_volume_source_energy_MeV(bad)
    bad,_,_ = _energy_observable_fixture(:ordinates;angular_count=2,weights=[0.25,1.75])
    bad.source.quadrature_weights[1] = Inf
    @test_throws ErrorException Radiant.get_volume_source_energy_MeV(bad)

    huge_values = fill(1.0e308,1,1,1)
    huge_norm = Source_Normalization(source_hash="overflow-test")
    huge_base = Anisotropic_Volume_Source(Photon(),[1],[1.0e300],[0.0,1.0e6],
        :isotropic,huge_values,huge_norm;source_tolerance=0.0)
    huge_wrapper = Energy_Moment_Volume_Source(huge_base,reshape(copy(huge_values),1,1,1,1))
    @test_throws ErrorException Radiant.get_volume_source_energy_MeV(huge_wrapper)
end
