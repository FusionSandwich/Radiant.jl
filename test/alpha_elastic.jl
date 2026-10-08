using Test
using Radiant

# These deliberately synthetic event DCS fixtures test transport mechanics.
# They do not qualify any ENDF evaluation or physical alpha scattering law.
function alpha_elastic_fixture(particle, Z, A; slope=0.0, identical=false,
                               spin=Radiant.is_alpha(particle) ? 0.0 : 0.5,
                               projectile_mass=Radiant.get_mass(particle),
                               projectile_charge=Int(Radiant.get_charge(particle)),target_charge=Z)
    E = 8.0
    mu = identical ? [0.0, 1.0] : [-1.0, 1.0]
    iso = Radiant.IsotopeElasticDB([E], mu, [[1.0+slope*x for x in mu]],
        1, identical, projectile_charge, target_charge, spin,
        projectile_mass, Radiant.elastic_target_mass(Z,A), Int[], Int[],
        nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    db = Radiant.ElasticScatteringENDFDB("database", "standard", true,
                                       Dict((Z,A) => iso))
    elastic = Elastic_Scattering()
    elastic.set_integrate_total_cs("numerical-gl")
    elastic.set_endf_db(Dict{Type,Radiant.ElasticScatteringENDFDB}(
        Radiant.get_type(particle) => db))
    return elastic, E / Radiant.M_E_C2_MEV
end

function alpha_elastic_feed(elastic, pin, pout, Z, A, Ei, process;
                            edges=[Ei,0.0], cutoff=Ei, order=2)
    return Radiant.feed_elastic_scattering([Z],[1.0],order,Ei,edges,length(edges)-1,
        elastic,1,1,Particle[pin,pout],process,pin,pout,[Ei,0.0],cutoff,
        false,false,[[A]],[[1.0]])
end

@testset "Alpha elastic kinematics: independent four-vector oracle" begin
    ma = Float64(Radiant.get_mass(Alpha()))
    mp = Float64(Radiant.get_mass(Proton()))
    @test Radiant.elastic_target_mass(1,1) == mp
    @test Radiant.elastic_target_mass(2,4) == ma
    # Swap projectile and target masses: an incorrect equal-mass shortcut fails.
    for (m,M) in ((ma,mp),(mp,ma),(ma,ma)), E in (0.05,8.0,200.0)
        s = m*m + M*M + 2*M*(m+E)
        root_s = sqrt(s)
        beta = sqrt(E*(E+2*m))/(m+E+M)
        gamma = (m+E+M)/root_s
        e1 = (s+m*m-M*M)/(2*root_s)
        e2 = (s+M*M-m*m)/(2*root_s)
        pstar = sqrt((s-(m+M)^2)*(s-(m-M)^2))/(2*root_s)
        for mu in (-1.0,-0.4,0.0,0.6,1.0)
            expected_s = gamma*(e1+beta*pstar*mu)-m
            expected_p = gamma*(e2-beta*pstar*mu)-M
            actual_s = Radiant.scattered_energy_relativistic(E,mu,m,M)
            actual_p = Radiant.recoil_energy_relativistic(E,mu,m,M)
            @test isapprox(actual_s,expected_s; atol=3e-10,rtol=1e-9)
            @test isapprox(actual_p,expected_p; atol=3e-10,rtol=1e-9)
            @test isapprox(actual_s+actual_p,E; rtol=2e-14)
            @test 0 <= actual_s <= E
            @test 0 <= actual_p <= E*(1+1e-14)
            @test isapprox(Radiant.mu_cm_from_Eout_relativistic(E,actual_s,m,M),mu; atol=2e-12)
            if -1 < mu < 1
                transverse = pstar*sqrt(1-mu*mu)
                z_s = gamma*(pstar*mu+beta*e1)
                z_p = gamma*(-pstar*mu+beta*e2)
                @test isapprox(Radiant.mu_lab_from_mu_cm(mu,E,m,M),z_s/hypot(z_s,transverse); atol=2e-11)
                @test isapprox(Radiant.recoil_mu_lab_from_mu_cm(mu,E,m,M),z_p/hypot(z_p,transverse); atol=2e-11)
            end
        end
    end
    @test Radiant.recoil_mu_lab_from_mu_cm(1.0,8.0,ma,mp) == 0.0
    @test_throws ArgumentError Radiant.scattered_energy_relativistic(0.0,0.0,ma,mp)
    @test_throws ArgumentError Radiant.scattered_energy_relativistic(-1.0,0.0,ma,mp)
    @test_throws ArgumentError Radiant.recoil_energy_relativistic(8.0,1.01,ma,mp)
    @test_throws ArgumentError Radiant.mu_cm_from_Eout_relativistic(8.0,NaN,ma,mp)
end

@testset "Alpha-H1 event number, angular direction and energy moments" begin
    pin, pout = Alpha(), Proton()
    elastic,Ei = alpha_elastic_fixture(pin,1,1; slope=0.5)
    Fs,Es = alpha_elastic_feed(elastic,pin,pin,1,1,Ei,"S")
    Fp,Ep = alpha_elastic_feed(elastic,pin,pout,1,1,Ei,"P")
    unit = Radiant.BARN_TO_CM2
    E = 8.0
    m, M = Radiant.get_mass(pin), Radiant.get_mass(pout)
    # Integrals on [-1,1] of q=1+mu/2: N=2 and int (1-mu)q=5/3.
    K = M*E*(E+2*m)/((m+M)^2+2*M*E)
    number = 4*pi*unit
    recoil_moment = (10*pi/3)*K/Radiant.M_E_C2_MEV*unit
    @test isapprox(sum(Fs[:,1]),number; rtol=2e-12)
    @test isapprox(sum(Fp[:,1]),number; rtol=2e-12)
    @test isapprox(sum(Ep),recoil_moment; rtol=2e-12)
    @test isapprox(sum(Es),Ei*number-recoil_moment; rtol=2e-12)
    @test isapprox(sum(Es)+sum(Ep),Ei*number; rtol=2e-12)
    @test isapprox(Radiant.tcs(elastic,Ei,Ei,pin,1; A=[1],atpercentA=[1.0]),number; rtol=2e-12)
    # Recoil energies above Ei/2 exist for alpha on hydrogen. The previous
    # Ei-Tmax upper bound and projectile energy inversion would lose this bin.
    Fhigh,_ = alpha_elastic_feed(elastic,pin,pout,1,1,Ei,"P"; edges=[Ei,0.5*Ei])
    @test Fhigh[1,1] > 0
    # Narrow event bin around mu=0: check the recoil lab angle independently.
    halfwidth = 1e-5
    edges = [K*(1+halfwidth),K*(1-halfwidth)]/Radiant.M_E_C2_MEV
    Fn,_ = alpha_elastic_feed(elastic,pin,pout,1,1,Ei,"P"; edges=edges)
    gamma = (E+m+M)/sqrt((m+M)^2+2*M*E)
    expected_mu = gamma/sqrt(1+gamma^2)
    @test isapprox(Fn[1,2]/Fn[1,1],expected_mu; atol=2e-10)
    @test abs(expected_mu-Radiant.mu_lab_from_mu_cm(0.0,E,m,M)) > 0.1
    # Energy moments must be invariant under an asymmetric redistribution of bins.
    _,Ebins = alpha_elastic_feed(elastic,pin,pout,1,1,Ei,"P"; edges=[Ei,0.57*Ei,0.19*Ei,0.0])
    @test isapprox(sum(Ebins),sum(Ep); rtol=2e-12)
    @test_throws ArgumentError alpha_elastic_feed(elastic,pin,pout,1,1,Ei,"S")
    Fwrong,_ = alpha_elastic_feed(elastic,pin,pin,1,1,Ei,"P")
    @test all(iszero,Fwrong) # H1 recoil is a proton, never an alpha.
    @test_throws ArgumentError alpha_elastic_feed(elastic,pin,pout,1,1,Ei,"P"; cutoff=-1.0)
    @test_throws ArgumentError alpha_elastic_feed(elastic,pin,pout,1,1,Ei,"P"; cutoff=Ei+1.0)
    @test_throws ErrorException alpha_elastic_feed(elastic,pin,pout,1,1,2*Ei,"P")
    @test_throws ArgumentError Radiant.tcs(elastic,0.0,0.0,pin,1; A=[1],atpercentA=[1.0])
end

@testset "Identical He4 inclusive S+P, spin and proton regression" begin
    for (pin,Z,A) in ((Alpha(),2,4),(Proton(),1,1))
        elastic,Ei = alpha_elastic_fixture(pin,Z,A; identical=true)
        Fs,Es = alpha_elastic_feed(elastic,pin,pin,Z,A,Ei,"S")
        Fp,Ep = alpha_elastic_feed(elastic,pin,pin,Z,A,Ei,"P")
        event_number = 2*pi*Radiant.BARN_TO_CM2
        @test isapprox(sum(Fs[:,1]),event_number; rtol=2e-12)
        @test isapprox(sum(Fp[:,1]),event_number; rtol=2e-12)
        @test isapprox(sum(Fs[:,1])+sum(Fp[:,1]),2*event_number; rtol=2e-12)
        @test isapprox(sum(Es),0.75*Ei*event_number; rtol=2e-12)
        @test isapprox(sum(Ep),0.25*Ei*event_number; rtol=2e-12)
        @test isapprox(sum(Es)+sum(Ep),Ei*event_number; rtol=2e-12)
        @test isapprox(Radiant.tcs(elastic,Ei,Ei,pin,Z; A=[A],atpercentA=[1.0]),event_number; rtol=2e-12)
        Sc,Sec = alpha_elastic_feed(elastic,pin,pin,Z,A,Ei,"S"; cutoff=0.8*Ei)
        Pc,Pec = alpha_elastic_feed(elastic,pin,pin,Z,A,Ei,"P"; cutoff=0.8*Ei)
        @test isapprox(sum(Sc[:,1]),0.6*event_number; rtol=2e-12)
        @test isapprox(sum(Pc[:,1]),sum(Sc[:,1]); rtol=2e-12)
        @test isapprox(sum(Sec)+sum(Pec),Ei*sum(Sc[:,1]); rtol=2e-12)
    end
    alpha = Alpha()
    elastic,Ei = alpha_elastic_fixture(alpha,2,4; identical=true)
    wrong,_ = alpha_elastic_feed(elastic,alpha,Proton(),2,4,Ei,"P")
    @test all(iszero,wrong) # He4 recoil is an alpha, never a proton.
    boson = Radiant.coulomb_dcs(0.3,1.0,0.2,0.01,true,0.0)
    fermion = Radiant.coulomb_dcs(0.3,1.0,0.2,0.01,true,0.5)
    @test boson > fermion > 0.0 # symmetry weight is distinct from two-output count.
    bad,Ei = alpha_elastic_fixture(alpha,2,4; identical=false)
    @test_throws ErrorException alpha_elastic_feed(bad,alpha,alpha,2,4,Ei,"S")
    bad,Ei = alpha_elastic_fixture(alpha,2,4; identical=true,spin=0.5)
    @test_throws ErrorException alpha_elastic_feed(bad,alpha,alpha,2,4,Ei,"S")
    bad,Ei = alpha_elastic_fixture(alpha,1,1; projectile_mass=Radiant.get_mass(Proton()))
    @test_throws ErrorException alpha_elastic_feed(bad,alpha,Proton(),1,1,Ei,"P")
    bad,Ei = alpha_elastic_fixture(alpha,1,1;projectile_charge=1)
    @test_throws ErrorException alpha_elastic_feed(bad,alpha,Proton(),1,1,Ei,"P")
    bad,Ei = alpha_elastic_fixture(alpha,1,1;target_charge=2)
    @test_throws ErrorException alpha_elastic_feed(bad,alpha,Proton(),1,1,Ei,"P")
end

@testset "Elastic configuration and missing data remain explicit" begin
    elastic = Elastic_Scattering()
    @test Radiant.get_types(elastic,Alpha,Alpha) == ["S","P"]
    @test Radiant.get_types(elastic,Alpha,Proton) == ["P"]
    elastic.set_interaction_types(Dict((Alpha,Proton)=>["P"]))
    @test Radiant.get_in_particles(elastic) == [Alpha]
    @test Radiant.get_out_particles(elastic) == [Proton]
    @test_throws ArgumentError elastic.set_interaction_types(Dict((Alpha,Proton)=>["S"]))
    @test_throws ArgumentError elastic.set_interaction_types(Dict((Alpha,Electron)=>["P"]))
    @test_throws KeyError alpha_elastic_feed(elastic,Alpha(),Proton(),1,1,8.0,"P")
    @test_throws ErrorException Radiant.init_elastic_scattering_endf_db(
        "__missing_alpha_fixture__",[(1,1)],Alpha(); data_root=joinpath(@__DIR__,"__absent__"))
end
