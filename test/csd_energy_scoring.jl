using Test
using Radiant

# Numerical scoring fixture, not a supplied-cross-section physical validation.
function manufactured_energy_score(method, dimension, full, energy_order, csd; imported=false)
    alpha,mat = Alpha("alpha"),Water()
    mat.set_density(1.25)
    cs = Cross_Sections()
    cs.set_particles([alpha]); cs.set_materials([mat]); cs.set_number_of_groups([1])
    mcs = Radiant.Multigroup_Cross_Sections(1)
    # S(x)=2+x on [-1,1], but an intentionally different stored average catches
    # failure to replace the group-averaged soft term in the legacy score.
    mcs.set_stopping_powers([7.0]); mcs.set_boundary_stopping_powers([3.0,1.0])
    mcs.set_energy_deposition([10.0,0.5]) # hard=3, soft=7; cutoff coefficient=.5
    if imported; mcs.stopping_powers=missing; end # FMAC-M/MATXS-style library
    cs.set_multigroup_cross_sections(reshape([mcs],1,1))
    geo=Geometry(); geo.set_dimension(dimension)
    for axis in ["x","y","z"][1:dimension]
        geo.set_number_of_regions(axis,1); geo.set_voxels_per_region(axis,[1])
        geo.set_region_boundaries(axis,[0.0,2.0])
        geo.set_boundary_conditions(uppercase(axis)*"-","void")
        geo.set_boundary_conditions(uppercase(axis)*"+","void")
    end
    geo.set_material_per_region(fill(mat,ntuple(_->1,dimension))); geo.build(cs)
    solver=method(); solver.set_particle(alpha); solver.set_solver_type(csd ? "BCSD" : "BTE")
    solver.set_is_full_coupling(full)
    for axis in ["x","y","z"][1:dimension]; solver.set_scheme(axis,"DG",2); end
    solver.set_scheme("E","DG",energy_order)
    solvers=Solvers(); solvers.add_solver(solver)
    sources=Fixed_Sources(cs,geo,solvers); sources.normalization_factor=4.0
    _,orders,nmom=solver.get_schemes(geo,full)
    f=zeros(1,1,nmom[5],1,1,1)
    f[1,1,1,1,1,1]=1.0
    if csd && energy_order>1
        f[1,1,2,1,1,1]=1/sqrt(3) # psi(x)=1+x
    end
    # Large zero-mean spatial moments must not enter the voxel average.
    for k in (orders[4]+1):nmom[5]; f[1,1,k,1,1,1]=100+k; end
    per=Radiant.Flux_Per_Particle(alpha); per.add_flux(f)
    if csd; per.add_flux_cutoff(fill(2.0,1,1,1,1,1)); end
    flux=Radiant.Flux(); flux.add_flux(per)
    return only(Radiant.energy_deposition(cs,geo,solvers,sources,flux,[alpha]))
end

@testset "CSD energy score contracts represented energy moments" begin
    # Independent physical-coordinate integral:
    # (1/2) integral[-1,1] (2+x)(1+x) dx = 7/3.
    for method in (SN,GN,DPN), dim in 1:3, full in (false,true), order in (1,2,3)
        soft = order==1 ? 2.0 : 7/3
        expected=(3.0+soft+0.5*2.0)/(1.25*4.0)
        @test manufactured_energy_score(method,dim,full,order,true) ≈ expected rtol=1e-14
    end
    for method in (SN,GN,DPN), dim in 1:3, full in (false,true)
        # Non-CSD keeps the original score; no stopping/cutoff correction applies.
        @test manufactured_energy_score(method,dim,full,1,false) == 10.0/(1.25*4.0)
        # Imported authoritative scores have no known separate Savg to replace.
        @test manufactured_energy_score(method,dim,full,2,true; imported=true) == 11.0/(1.25*4.0)
    end
end
