using Radiant
using SHA
using DelimitedFiles
using Printf

length(ARGS) in (4,5,6,7,8) || error("usage: julia run_radiant.jl <stopping.csv> <geant4_scores.csv> <voxels-per-material> <lower-energy-groups> [energy-MeV [slab-cm [output.csv [lower-energy-MeV]]]]")
table_path, g4_scores_path = ARGS[1:2]
nvox = parse(Int,ARGS[3])
ng = parse(Int,ARGS[4])
energy_filter = length(ARGS) == 5 ? parse(Float64,ARGS[5]) : nothing
if length(ARGS) >= 6
    energy_filter = parse(Float64,ARGS[5])
end
nvox > 0 && ng > 0 || error("Voxel and energy-group counts must be positive.")

table = readdlm(table_path,',',Float64;skipstart=1)
g4 = readdlm(g4_scores_path,',',Float64;skipstart=1)
energies = vec(table[:,1])
all(diff(energies) .> 0) || error("Geant4 stopping grid is not increasing.")
all(isfinite,table) && all(table[:,2:3] .> 0) || error("Geant4 stopping table has invalid values.")
table_hash = bytes2hex(sha256(read(table_path)))
slab_cm = length(ARGS) >= 6 ? parse(Float64,ARGS[6]) : 0.01
slab_cm > 0 || error("Slab thickness must be positive.")
out_path = length(ARGS) >= 7 ? ARGS[7] : joinpath(dirname(table_path),"radiant_scores_mass_corrected.csv")
lower_override = length(ARGS) == 8 ? parse(Float64,ARGS[8]) : nothing
legacy_header = "energy_MeV,voxels_per_material,lower_energy_groups,Al_edep_MeV,Cu_edep_MeV,source_rate,solve_s,classification,table_sha256"
header_with_slab = legacy_header * ",slab_cm"
header_with_lower = header_with_slab * ",lower_MeV"
if !isfile(out_path)
    open(out_path,"w") do io
        println(io,header_with_lower)
    end
end
out_header = open(io -> readline(io),out_path)
out_header in (legacy_header,header_with_slab,header_with_lower) ||
    error("Unknown Radiant score-file schema; refusing to append.")

for row in eachrow(g4)
    E = row[1]
    if !isnothing(energy_filter) && !isapprox(E,energy_filter;atol=1e-9)
        continue
    end
    densities = (row[3],row[4])
    lower = isnothing(lower_override) ? (E <= 20 ? 1.0 : E <= 200 ? E-5.0 : E-2.0) : lower_override
    1.0 <= lower < E-0.005 || error("Invalid Radiant lower-energy boundary.")
    boundaries = vcat([E+0.005],collect(range(E-0.005,stop=lower,length=ng+1)))
    proton = Proton()
    materials = Material[]
    material_data = Proton_Material_Data[]
    nonelastic = Proton_Nonelastic_Data[]
    for (i,id) in enumerate(("Al","Cu"))
        material = Material(id)
        Radiant.set_density(material,densities[i])
        model = Tabulated_Ion_Transport_Model(
            species_id="proton",material_id=id,
            energy_MeV=energies,
            electronic_stopping_MeV_cm=vec(table[:,i+1]),
            nuclear_stopping_MeV_cm=zeros(length(energies)),
            energy_straggling_variance_MeV2_cm=zeros(length(energies)),
            angular_variance_rad2_cm=zeros(length(energies)),
            data_hash=table_hash,qualification_status=:candidate)
        push!(materials,material)
        push!(material_data,Proton_Material_Data(material,model;
            material_state="solid",density_g_cm3=densities[i],
            source_sha256=table_hash,source_path=table_path,
            uncertainty_status=:unquantified))
        nhash = bytes2hex(sha256("$id/no-nonelastic-in-controlled-G4-model"))
        push!(nonelastic,Proton_Nonelastic_Data(
            material_id=id,energy_MeV=[lower,E+0.005],
            removal_cm_inv=[0.0,0.0],event_family_ids=String[],
            source_sha256=nhash,qualification_status=:synthetic))
    end
    binding = bind_proton_native(proton,material_data,nonelastic,boundaries;
        qualification_mode=:software)
    cs = binding.cross_sections
    Radiant.build(cs)
    geometry = Geometry()
    Radiant.set_dimension(geometry,1)
    Radiant.set_number_of_regions(geometry,"x",2)
    Radiant.set_region_boundaries(geometry,"x",[0.0,slab_cm,2*slab_cm])
    Radiant.set_voxels_per_region(geometry,"x",[nvox,nvox])
    Radiant.set_boundary_conditions(geometry,"x-","void")
    Radiant.set_boundary_conditions(geometry,"x+","void")
    Radiant.set_material_per_region(geometry,materials)
    Radiant.build(geometry,cs)
    solver = SN()
    Radiant.set_particle(solver,proton)
    Radiant.set_solver_type(solver,"CSD")
    Radiant.set_quadrature(solver,"gauss-legendre",2)
    Radiant.set_legendre_order(solver,0)
    Radiant.set_angular_boltzmann(solver,"galerkin-d")
    Radiant.set_scheme(solver,"x","DD",1)
    Radiant.set_scheme(solver,"E","DG",1)
    solvers = Solvers()
    Radiant.add_solver(solvers,solver)
    normalization = Source_Normalization(
        basis=:per_source_particle,source_rate_per_s=1.0,
        source_hash=bytes2hex(sha256("controlled-g4-source-E=$E")),
        provenance=Dict("case"=>"geant4-ionisation-only-slab"))
    source_values = zeros(nvox,length(boundaries)-1,1)
    source_values[:,end,1] .= 1/slab_cm
    source = Anisotropic_Volume_Source(proton,collect(1:nvox),
        fill(slab_cm/nvox,nvox),reverse(boundaries).*1.0e6,
        :isotropic,source_values,normalization)
    sources = Fixed_Sources(cs,geometry,solvers)
    Radiant.add_source(sources,source)
    cu = Computation_Unit()
    Radiant.set_cross_sections(cu,cs)
    Radiant.set_geometry(cu,geometry)
    Radiant.set_solvers(cu,solvers)
    Radiant.set_sources(cu,sources)
    field = Electromagnetic_Field()
    Radiant.set_magnetic_field(field,[0.0,0.0,0.0])
    Radiant.set_electromagnetic_field(cu,field)
    source_rate = sum(get_volume_source_rate(source))
    isapprox(source_rate,1.0;rtol=1e-12) || error("Source normalization is not one proton per source.")
    seconds = @elapsed Radiant.run(cu)
    score = vec(Radiant.get_energy_deposition(cu,proton))
    length(score) == 2*nvox || error("Unexpected Radiant deposition shape.")
    # Radiant returns a per-voxel mass-normalized score, not deposited MeV.
    # In 1D its documented units include a cm factor; multiply by the
    # material density and voxel width for MeV per source proton per unit area.
    al = sum(score[1:nvox])*densities[1]*(slab_cm/nvox)
    copper = sum(score[nvox+1:end])*densities[2]*(slab_cm/nvox)
    classification = proton_binding_receipt(binding)["classification"]
    @printf("E=%.5f MeV nvox=%d ng=%d Al=%.9g Cu=%.9g solve=%.3f s %s\n",E,nvox,ng,al,copper,seconds,classification)
    open(out_path,"a") do io
        @printf(io,"%.8f,%d,%d,%.12g,%.12g,%.12g,%.6f,%s,%s",E,nvox,ng,al,copper,source_rate,seconds,classification,table_hash)
        out_header != legacy_header && @printf(io,",%.9g",slab_cm)
        out_header == header_with_lower && @printf(io,",%.9g",lower)
        println(io)
    end
end
