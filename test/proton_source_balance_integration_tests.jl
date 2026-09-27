using Test, LinearAlgebra

# public_energy_problem and public_energy_base are the registered public API
# fixtures. Two cells and spatial DG2 challenge the diagnostic's extraction of
# spatial averages and cancellation of the unobserved internal spatial face.
@testset "Public volume source and native spatial DG2 weak balance" begin
    nodes = [-sqrt(3/5),0.0,sqrt(3/5)]
    weights = [5/9,8/9,5/9]
    edges = [1e6,5e6,10e6]
    for oe in 1:4, full in (false,true)
        problem = public_energy_problem(;energy_order=oe,full=full)
        powers = zeros(2,2,1,oe)
        for v in 1:2,g in 1:2
            powers[v,g,1,:] .= (v+g)*1e-6 .* [1.0,-0.2,0.3,0.1][1:oe]
        end
        moments = source_energy_moments(edges,powers;energy_order=oe)
        source = Energy_Moment_Volume_Source(public_energy_base(problem,moments),moments)
        sources = Fixed_Sources(problem.cs,problem.geo,problem.solvers)
        Radiant.add_source(sources,source)
        Radiant.build(sources)
        flux = Radiant.transport(problem.cs,problem.geo,problem.solvers,sources;
                                 retain_boundary_flux=true)

        # Independent Gauss3 integral is exact for E times degree <=3 source.
        injection = 0.0
        for v in 1:2,g in 1:2
            width = edges[g+1]-edges[g]
            mid = (edges[g+1]+edges[g])/2
            injection += source.source.voxel_volumes_cm3[v]*width/2*sum(
                weights[d]*(mid+width/2*nodes[d])*
                sum(powers[v,g,1,k]*nodes[d]^(k-1) for k in 1:oe)
                for d in eachindex(nodes))/1e6
        end
        @test sum(get_volume_source_energy_MeV(source)) ≈ injection atol=1e-10 rtol=1e-13
        diagnostic = proton_discrete_energy_balance(
            problem.binding,problem.geo,problem.solvers,sources,flux;
            injected_energy_MeV=injection)
        @test diagnostic.orders == [2,1,1,oe]
        @test diagnostic.fully_coupled == full
        @test diagnostic.convergence_verified
        @test diagnostic.energy_test_represented == (oe >= 2)
        @test abs(diagnostic.injection_projection_difference_MeV) < 1e-9
        @test maximum(abs,diagnostic.particle_balance_residual_by_group) < 1e-9
        @test maximum(abs,diagnostic.tested_balance_residual_by_group_MeV) < 1e-9
        @test maximum(abs,diagnostic.stored_minus_reconstructed_cutoff) < 1e-9
        @test abs(diagnostic.weak_balance_residual_MeV) < 1e-9
        @test !diagnostic.strict_ledger_accepted && !diagnostic.response_score_replaced
        @test !diagnostic.physical_validation && !diagnostic.publication_ready && !diagnostic.pr_ready
        println("PUBLIC_SOURCE_SPATIAL_DG2_BALANCE energy_order=",oe," full=",full,
                " independent_injection_MeV=",injection,
                " weak_residual_MeV=",diagnostic.weak_balance_residual_MeV,
                " midpoint_response_residual_MeV=",diagnostic.represented_midpoint_response_residual_MeV,
                " minimum_scalar_cutoff=",diagnostic.minimum_stored_scalar_cutoff,
                " physical_validation=false")
    end
end
