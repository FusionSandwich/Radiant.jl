using Radiant
using LinearAlgebra
using Test

const SOFTWARE_CLASSIFICATION = "SOFTWARE_VERIFICATION_XS"
const ENERGY_BOUNDS_ONE_GROUP_MEV = [1.0,0.001]

function make_material(tag::String,density_g_cm3::Float64=1.0)
    material = Material(tag)
    Radiant.set_density(material,density_g_cm3)
    return material
end

function make_uniform_em_library(
    materials::Vector{Material};
    number_of_groups::Int=1,
    include_electron::Bool=false,
    photon_total::Vector{Float64},
    photon_absorption::Vector{Float64},
    photon_deposition::Vector{Float64},
    photon_to_electron::Vector{Float64}=zeros(length(materials)),
    electron_total::Vector{Float64}=fill(0.7,length(materials)),
    electron_absorption::Vector{Float64}=fill(0.7,length(materials)),
    electron_deposition::Vector{Float64}=fill(0.45,length(materials)),
    electron_charge_deposition::Vector{Float64}=fill(0.20,length(materials)),
)
    nmaterial = length(materials)
    all(length(values) == nmaterial for values in (
        photon_total,photon_absorption,photon_deposition,photon_to_electron,
        electron_total,electron_absorption,electron_deposition,electron_charge_deposition,
    )) || error("Every material-response vector must match the material count.")

    photon = Photon()
    particles = include_electron ? Particle[photon,Electron()] : Particle[photon]
    nparticle = length(particles)
    ngroup = Int(number_of_groups)
    ngroup ≥ 1 || error("At least one energy group is required.")

    # High-to-low Radiant convention. The exact values are synthetic but stable and positive.
    energy_boundaries = ngroup == 1 ? copy(ENERGY_BOUNDS_ONE_GROUP_MEV) :
        collect(10.0 .^ range(0.0,-3.0;length=ngroup+1))

    libraries = Array{Radiant.Multigroup_Cross_Sections}(undef,nparticle,nmaterial)
    for material_index in 1:nmaterial
        gamma = Radiant.Multigroup_Cross_Sections(ngroup)
        Radiant.set_total(gamma,fill(photon_total[material_index],ngroup))
        Radiant.set_absorption(gamma,fill(photon_absorption[material_index],ngroup))
        Radiant.set_boundary_stopping_powers(gamma,zeros(ngroup+1))
        Radiant.set_stopping_powers(gamma,zeros(ngroup))
        Radiant.set_momentum_transfer(gamma,zeros(ngroup))
        Radiant.set_energy_deposition(
            gamma,vcat(fill(photon_deposition[material_index],ngroup),0.0),
        )
        Radiant.set_charge_deposition(gamma,zeros(ngroup+1))

        gamma_to_gamma = zeros(Float64,ngroup,ngroup,1)
        Radiant.set_scattering(gamma,gamma_to_gamma)
        if include_electron
            gamma_to_electron = zeros(Float64,ngroup,ngroup,1)
            for group in 1:ngroup
                gamma_to_electron[group,group,1] = photon_to_electron[material_index]
            end
            Radiant.set_scattering(gamma,gamma_to_electron)
        end
        libraries[1,material_index] = gamma

        if include_electron
            electron = Radiant.Multigroup_Cross_Sections(ngroup)
            Radiant.set_total(electron,fill(electron_total[material_index],ngroup))
            Radiant.set_absorption(electron,fill(electron_absorption[material_index],ngroup))
            Radiant.set_boundary_stopping_powers(electron,zeros(ngroup+1))
            Radiant.set_stopping_powers(electron,zeros(ngroup))
            Radiant.set_momentum_transfer(electron,zeros(ngroup))
            Radiant.set_energy_deposition(
                electron,vcat(fill(electron_deposition[material_index],ngroup),0.0),
            )
            Radiant.set_charge_deposition(
                electron,vcat(fill(electron_charge_deposition[material_index],ngroup),0.0),
            )
            Radiant.set_scattering(electron,zeros(Float64,ngroup,ngroup,1)) # e -> gamma
            Radiant.set_scattering(electron,zeros(Float64,ngroup,ngroup,1)) # e -> e
            libraries[2,material_index] = electron
        end
    end

    cross_sections = Cross_Sections()
    Radiant.set_materials(cross_sections,materials)
    Radiant.set_particles(cross_sections,particles)
    Radiant.set_number_of_groups(cross_sections,fill(ngroup,nparticle))
    Radiant.set_energy_boundaries(
        cross_sections,[copy(energy_boundaries) for _ in 1:nparticle],
    )
    Radiant.set_multigroup_cross_sections(cross_sections,libraries)
    cross_sections.is_build = true
    return cross_sections,particles,energy_boundaries
end

function make_geometry(
    cross_sections::Cross_Sections,
    materials::Vector{Material},
    widths_cm::Vector{Float64},
    voxels_per_region::Vector{Int64},
)
    length(materials) == length(widths_cm) == length(voxels_per_region) || error(
        "Materials, widths, and voxel counts must have the same length.",
    )
    boundaries = zeros(Float64,length(widths_cm)+1)
    for index in eachindex(widths_cm)
        widths_cm[index] > 0.0 || error("Region widths must be positive.")
        boundaries[index+1] = boundaries[index] + widths_cm[index]
    end

    geometry = Geometry()
    Radiant.set_type(geometry,"cartesian")
    Radiant.set_dimension(geometry,1)
    Radiant.set_number_of_regions(geometry,"x",length(materials))
    Radiant.set_region_boundaries(geometry,"x",boundaries)
    Radiant.set_voxels_per_region(geometry,"x",voxels_per_region)
    Radiant.set_boundary_conditions(geometry,"x-","void")
    Radiant.set_boundary_conditions(geometry,"x+","void")
    Radiant.set_material_per_region(geometry,materials)
    Radiant.build(geometry,cross_sections)
    return geometry
end

function make_solvers(
    particles::Vector{Particle},
    quadrature_order::Int;
    maximum_generations::Int=8,
)
    solvers = Solvers()
    for particle in particles
        solver = SN()
        Radiant.set_particle(solver,particle)
        Radiant.set_solver_type(solver,"BTE")
        Radiant.set_quadrature(solver,"gauss-legendre",quadrature_order)
        Radiant.set_legendre_order(solver,0)
        Radiant.set_angular_boltzmann(solver,"galerkin-d")
        Radiant.set_convergence_criterion(solver,1.0e-11)
        Radiant.set_maximum_iteration(solver,500)
        Radiant.set_scheme(solver,"x","DD",1)
        Radiant.add_solver(solvers,solver)
    end
    Radiant.set_maximum_number_of_generations(solvers,maximum_generations)
    Radiant.set_convergence_criterion(solvers,1.0e-10)
    Radiant.set_convergence_type(solvers,"flux")
    return solvers
end

function make_boundary_source(
    particle::Particle,
    cross_sections::Cross_Sections,
    geometry::Geometry,
    solver::SN,
    energy_boundaries_mev::Vector{Float64};
    angular_mode::Symbol=:isotropic_incoming,
    split_group_current::Bool=true,
    source_hash::String="synthetic-opensn-source",
)
    _,weights,directions,_ = Radiant._solver_quadrature(solver,geometry)
    normal = [-1.0 0.0 0.0]
    incoming = [
        direction for direction in axes(directions,1)
        if dot(view(directions,direction,:),view(normal,1,:)) < 0.0
    ]
    isempty(incoming) && error("The selected quadrature has no X- incoming directions.")

    ngroup = Radiant.get_number_of_groups(cross_sections,particle)
    directional_current = zeros(Float64,1,ngroup,length(weights))
    group_scale = split_group_current ? 1.0/ngroup : 1.0

    if angular_mode == :isotropic_incoming
        factors = [
            weights[direction] * abs(dot(
                view(directions,direction,:),view(normal,1,:),
            )) for direction in incoming
        ]
        denominator = sum(factors)
        denominator > 0.0 || error("Incoming-current quadrature factor is zero.")
        for group in 1:ngroup, (local_index,direction) in enumerate(incoming)
            directional_current[1,group,direction] = group_scale * factors[local_index]/denominator
        end
    elseif angular_mode == :near_normal
        direction = incoming[argmax([directions[index,1] for index in incoming])]
        directional_current[1,:,direction] .= group_scale
    elseif angular_mode == :oblique
        positive = sort(incoming;by=index -> directions[index,1])
        direction = positive[max(1,cld(length(positive),2))]
        directional_current[1,:,direction] .= group_scale
    else
        error("Unknown angular source mode $(angular_mode).")
    end

    normalization = Source_Normalization(
        basis=:per_history,
        source_rate_per_s=1.0,
        symmetry_factor=1.0,
        time_class=:prompt,
        source_hash=source_hash,
        provenance=Dict(
            "producer" => "synthetic-opensn-replay",
            "classification" => SOFTWARE_CLASSIFICATION,
        ),
    )
    source_edges_eV = reverse(1.0e6 .* energy_boundaries_mev)
    return boundary_source_from_directional_current(
        particle,
        [1],
        [0.0 0.0 0.0],
        [1.0],
        normal,
        [0.0 1.0 0.0],
        [0.0 0.0 -1.0],
        source_edges_eV,
        directions,
        weights,
        directional_current,
        normalization;
        provenance=Dict(
            "source_contract" => "synthetic-opensn-angular-current/v1",
            "classification" => SOFTWARE_CLASSIFICATION,
        ),
    )
end

function execute_case(;
    case_id::String,
    materials::Vector{Material},
    widths_cm::Vector{Float64},
    voxels_per_region::Vector{Int64},
    number_of_groups::Int=1,
    quadrature_order::Int=8,
    angular_mode::Symbol=:isotropic_incoming,
    include_electron::Bool=false,
    photon_total::Vector{Float64},
    photon_absorption::Vector{Float64},
    photon_deposition::Vector{Float64},
    photon_to_electron::Vector{Float64}=zeros(length(materials)),
    electron_total::Vector{Float64}=fill(0.7,length(materials)),
    electron_absorption::Vector{Float64}=fill(0.7,length(materials)),
    electron_deposition::Vector{Float64}=fill(0.45,length(materials)),
    electron_charge_deposition::Vector{Float64}=fill(0.20,length(materials)),
)
    cross_sections,particles,energy_boundaries = make_uniform_em_library(
        materials;
        number_of_groups=number_of_groups,
        include_electron=include_electron,
        photon_total=photon_total,
        photon_absorption=photon_absorption,
        photon_deposition=photon_deposition,
        photon_to_electron=photon_to_electron,
        electron_total=electron_total,
        electron_absorption=electron_absorption,
        electron_deposition=electron_deposition,
        electron_charge_deposition=electron_charge_deposition,
    )
    geometry = make_geometry(cross_sections,materials,widths_cm,voxels_per_region)
    solvers = make_solvers(particles,quadrature_order)
    photon_solver = Radiant.get_method(solvers,particles[1])
    source = make_boundary_source(
        particles[1],cross_sections,geometry,photon_solver,energy_boundaries;
        angular_mode=angular_mode,
        source_hash="$(case_id)-source",
    )
    fixed_sources = Fixed_Sources(cross_sections,geometry,solvers)
    Radiant.add_source(fixed_sources,source)

    calculation = Computation_Unit()
    Radiant.set_cross_sections(calculation,cross_sections)
    Radiant.set_geometry(calculation,geometry)
    Radiant.set_solvers(calculation,solvers)
    Radiant.set_sources(calculation,fixed_sources)
    Radiant.run(calculation)

    photon_flux = vec(Radiant.get_flux(calculation,particles[1]))
    photon_deposition_result = vec(Radiant.get_energy_deposition(calculation,particles[1]))
    all(isfinite,photon_flux) || error("$(case_id) produced nonfinite photon flux.")
    all(value -> value ≥ -1.0e-12,photon_flux) || error("$(case_id) produced negative photon flux.")
    all(isfinite,photon_deposition_result) || error("$(case_id) produced nonfinite photon deposition.")

    electron_flux_sum = 0.0
    electron_deposition_sum = 0.0
    charge_deposition_sum = 0.0
    if include_electron
        electron_flux = vec(Radiant.get_flux(calculation,particles[2]))
        electron_deposition_result = vec(Radiant.get_energy_deposition(calculation,particles[2]))
        charge_deposition_result = vec(Radiant.get_charge_deposition(calculation,particles[2]))
        all(isfinite,electron_flux) || error("$(case_id) produced nonfinite electron flux.")
        all(value -> value ≥ -1.0e-12,electron_flux) || error("$(case_id) produced negative electron flux.")
        electron_flux_sum = sum(electron_flux)
        electron_deposition_sum = sum(electron_deposition_result)
        charge_deposition_sum = sum(charge_deposition_result)
    end

    receipt = get_projection_receipts(fixed_sources)[1]
    return (
        case_id=case_id,
        photon_flux=photon_flux,
        photon_flux_sum=sum(photon_flux),
        photon_deposition=photon_deposition_result,
        photon_deposition_sum=sum(photon_deposition_result),
        electron_flux_sum=electron_flux_sum,
        electron_deposition_sum=electron_deposition_sum,
        charge_deposition_sum=charge_deposition_sum,
        source_current=sum(receipt.target_current),
        projected_current=sum(receipt.projected_current),
        current_error=receipt.max_relative_error,
        geometry=geometry,
        source=source,
        calculation=calculation,
    )
end

function relative_difference(first::Real,second::Real)
    return abs(Float64(first)-Float64(second))/max(abs(Float64(second)),1.0e-14)
end

@testset "Radiant executed EM qualification" begin
    vacuum_material = [make_material("vacuum-like")]
    vacuum = execute_case(
        case_id="photon-vacuum",
        materials=vacuum_material,
        widths_cm=[1.0],
        voxels_per_region=[32],
        photon_total=[0.0],
        photon_absorption=[0.0],
        photon_deposition=[0.0],
    )
    @test maximum(vacuum.photon_flux)-minimum(vacuum.photon_flux) ≤ 1.0e-9
    @test vacuum.photon_deposition_sum == 0.0
    @test vacuum.current_error ≤ 1.0e-10

    absorber_material = [make_material("absorber")]
    absorber = execute_case(
        case_id="single-material-photon",
        materials=absorber_material,
        widths_cm=[1.0],
        voxels_per_region=[64],
        angular_mode=:near_normal,
        photon_total=[0.8],
        photon_absorption=[0.8],
        photon_deposition=[0.55],
    )
    @test absorber.photon_flux[1] > absorber.photon_flux[end]
    @test absorber.photon_deposition_sum > 0.0

    two_materials = [make_material("front"),make_material("rear")]
    two_layer = execute_case(
        case_id="two-material-photon",
        materials=two_materials,
        widths_cm=[0.4,0.6],
        voxels_per_region=[24,36],
        photon_total=[0.25,1.0],
        photon_absorption=[0.25,1.0],
        photon_deposition=[0.18,0.70],
    )
    @test two_layer.photon_deposition_sum > 0.0
    @test two_layer.photon_flux[1] > two_layer.photon_flux[end]

    layer_names = ["cu-front","ag","ybco","buffer","hastelloy","cu-rear","solder","insulation"]
    stack_materials = [make_material(name) for name in layer_names]
    stack_widths = [0.002,0.0002,0.0001,0.00002,0.005,0.002,0.002,0.005]
    stack_xs = [25.0,40.0,80.0,60.0,20.0,25.0,30.0,10.0]
    stack_deposition = 0.65 .* stack_xs
    base_voxels = Int64[4,2,4,2,8,4,4,8]

    stack_normal = execute_case(
        case_id="eight-layer-near-normal",
        materials=stack_materials,
        widths_cm=stack_widths,
        voxels_per_region=base_voxels,
        angular_mode=:near_normal,
        photon_total=stack_xs,
        photon_absorption=stack_xs,
        photon_deposition=stack_deposition,
    )
    stack_oblique = execute_case(
        case_id="eight-layer-oblique",
        materials=stack_materials,
        widths_cm=stack_widths,
        voxels_per_region=base_voxels,
        angular_mode=:oblique,
        photon_total=stack_xs,
        photon_absorption=stack_xs,
        photon_deposition=stack_deposition,
    )
    @test stack_normal.photon_deposition_sum > 0.0
    @test stack_oblique.photon_deposition_sum > stack_normal.photon_deposition_sum
    ybco_indices = (sum(base_voxels[1:2])+1):sum(base_voxels[1:3])
    @test sum(stack_normal.photon_deposition[ybco_indices]) > 0.0

    mixed = execute_case(
        case_id="mixed-photon-electron",
        materials=[make_material("em-medium")],
        widths_cm=[1.0],
        voxels_per_region=[64],
        include_electron=true,
        angular_mode=:near_normal,
        photon_total=[0.9],
        photon_absorption=[0.65],
        photon_deposition=[0.30],
        photon_to_electron=[0.25],
        electron_total=[1.2],
        electron_absorption=[1.2],
        electron_deposition=[0.75],
        electron_charge_deposition=[0.35],
    )
    @test mixed.photon_deposition_sum > 0.0
    @test mixed.electron_flux_sum > 0.0
    @test mixed.electron_deposition_sum > 0.0
    @test mixed.charge_deposition_sum > 0.0

    spatial = [
        execute_case(
            case_id="spatial-$(factor)",
            materials=stack_materials,
            widths_cm=stack_widths,
            voxels_per_region=Int64.(factor .* base_voxels),
            angular_mode=:isotropic_incoming,
            photon_total=stack_xs,
            photon_absorption=stack_xs,
            photon_deposition=stack_deposition,
        ) for factor in (1,2,4)
    ]
    @test relative_difference(
        spatial[2].photon_deposition_sum,spatial[3].photon_deposition_sum,
    ) < relative_difference(
        spatial[1].photon_deposition_sum,spatial[2].photon_deposition_sum,
    )
    @test relative_difference(spatial[2].photon_deposition_sum,spatial[3].photon_deposition_sum) < 0.03

    angular = [
        execute_case(
            case_id="angular-$(order)",
            materials=stack_materials,
            widths_cm=stack_widths,
            voxels_per_region=Int64.(2 .* base_voxels),
            quadrature_order=order,
            angular_mode=:isotropic_incoming,
            photon_total=stack_xs,
            photon_absorption=stack_xs,
            photon_deposition=stack_deposition,
        ) for order in (4,8,16)
    ]
    @test relative_difference(
        angular[2].photon_deposition_sum,angular[3].photon_deposition_sum,
    ) < relative_difference(
        angular[1].photon_deposition_sum,angular[2].photon_deposition_sum,
    )
    @test relative_difference(angular[2].photon_deposition_sum,angular[3].photon_deposition_sum) < 0.03

    energy_one = execute_case(
        case_id="energy-one-group",
        materials=absorber_material,
        widths_cm=[1.0],
        voxels_per_region=[64],
        number_of_groups=1,
        photon_total=[0.8],
        photon_absorption=[0.8],
        photon_deposition=[0.55],
    )
    energy_two = execute_case(
        case_id="energy-two-group-split",
        materials=absorber_material,
        widths_cm=[1.0],
        voxels_per_region=[64],
        number_of_groups=2,
        photon_total=[0.8],
        photon_absorption=[0.8],
        photon_deposition=[0.55],
    )
    @test relative_difference(
        energy_one.photon_deposition_sum,energy_two.photon_deposition_sum,
    ) < 1.0e-8

    println("RADIANT_EM_PASS|classification=$(SOFTWARE_CLASSIFICATION)")
    println("RADIANT_EM_PASS|mixed_electron_flux=$(mixed.electron_flux_sum)")
    println("RADIANT_EM_PASS|mixed_electron_deposition=$(mixed.electron_deposition_sum)")
    println("RADIANT_EM_PASS|mixed_charge_deposition=$(mixed.charge_deposition_sum)")
    println("RADIANT_EM_PASS|oblique_to_normal_deposition_ratio=$(stack_oblique.photon_deposition_sum/stack_normal.photon_deposition_sum)")
    println("OPENSN_RADIANT_SURFACE_REPLAY_PASS|max_current_error=$(maximum(case.current_error for case in (vacuum,absorber,two_layer,stack_normal,stack_oblique,mixed)))")
end
