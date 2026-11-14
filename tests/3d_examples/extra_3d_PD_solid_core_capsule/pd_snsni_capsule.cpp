/** Frozen S-NSNI/RKPM DPC core + existing layered NOSB-PD energy-damage shell.
 * One common explicit displacement/velocity state, physical lumped masses,
 * common-trace vector penalty and frictionless reconstructed outer contact.
 * No WCSPH fluid, material retuning, particle erosion or force rescaling.
 */
#include "pd_snsni_capsule_fixture.h"
#include "solid_snsni_plastic.h"
#include "pd_radial_layered_surface.h"
#include "recovery_checkpoint.h"
using namespace SPH::solid_dynamics;
using namespace SPH::peridynamics;
using Field=StdVec<Vecd>;

class SnsniCoupledSurface:public ReferenceSolidSurface
{
 public:
    SnsniCoupledSurface(SolidForcePort&p,const StdVec<Vecd>&points,const SnsniReference&r)
      :ReferenceSolidSurface(p,points,1,DeferredBuildTag{})
    {
        for(const auto&q:points){Row row;Real sum=0;Vecd affine=Vecd::Zero();
            for(const auto&e:r.shape(q)){row.index.push_back(e.j);row.weight.push_back(e.w);sum+=e.w;affine+=e.w*r.X[e.j];}
            affine_error=std::max({affine_error,std::abs(sum-1),(affine-q).norm()/r.dp});rows.push_back(std::move(row));}
        if(affine_error>1e-8)throw std::runtime_error("Core trace failed affine reproduction");
    }
};

void vtkPoints(const fs::path&path,const Field&x,const Field&v,const StdVec<Real>&damage,
               const StdVec<Real>&pressure,const StdVec<Real>&J,const StdVec<Real>&plastic,
               const StdVec<Real>&dissipation,const Field&normal)
{
    std::ofstream f(path);f<<std::setprecision(17);size_t n=x.size();
    f<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\"><PolyData><Piece NumberOfPoints=\""<<n<<"\" NumberOfVerts=\""<<n<<"\"><Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(const auto&a:x)f<<a[0]<<' '<<a[1]<<' '<<a[2]<<'\n';
    f<<"</DataArray></Points><Verts><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">";for(size_t i=0;i<n;++i)f<<i<<' ';
    f<<"</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">";for(size_t i=0;i<n;++i)f<<i+1<<' ';f<<"</DataArray></Verts><PointData>";
    auto scalar=[&](const char*name,const StdVec<Real>&a){if(a.empty())return;f<<"<DataArray type=\"Float64\" Name=\""<<name<<"\" format=\"ascii\">";for(Real z:a)f<<z<<' ';f<<"</DataArray>";};
    auto vector=[&](const char*name,const Field&a){if(a.empty())return;f<<"<DataArray type=\"Float64\" Name=\""<<name<<"\" NumberOfComponents=\"3\" format=\"ascii\">";for(const auto&z:a)f<<z[0]<<' '<<z[1]<<' '<<z[2]<<' ';f<<"</DataArray>";};
    scalar("Damage",damage);scalar("Pressure_Pa",pressure);scalar("J",J);scalar("PlasticShear",plastic);scalar("FractureDissipationDensity",dissipation);vector("Velocity",v);vector("CrackNormalReference",normal);
    f<<"</PointData></Piece></PolyData></VTKFile>\n";if(!f)throw std::runtime_error("VTP output failed");
}

void vtkPlates(const fs::path&path,Real R,Real compression)
{
    std::ofstream f(path);f<<std::setprecision(17)<<"<VTKFile type=\"PolyData\" version=\"0.1\"><PolyData><Piece NumberOfPoints=\"8\" NumberOfPolys=\"2\"><Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
    for(int side:{-1,1})for(auto xy:{std::pair{-1,-1},std::pair{1,-1},std::pair{1,1},std::pair{-1,1}})f<<xy.first*1.35*R<<' '<<xy.second*1.35*R<<' '<<side*(R-compression/2)<<' ';
    f<<"</DataArray></Points><Polys><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">0 1 2 3 4 5 6 7</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">4 8</DataArray></Polys></Piece></PolyData></VTKFile>\n";
}

int runSnsniCapsule(const HybridOptions&o,const std::string&cells)
{
    if(!o.layered_kinematics||o.shell_material!="finite_volume_biot"||o.fracture_model!="crack_band_energy"||o.plate_friction!=0)
        throw std::invalid_argument("This entry requires layered NOSB, finite-volume Biot crack-band damage and frictionless contact");
    fs::path out=fs::absolute(o.output);if(fs::exists(out)&&!fs::is_empty(out))throw std::runtime_error("Output must be new or empty; no retry/overwrite");
    fs::create_directories(out/"vtp");const Real R=o.diameter/2,rc=R-o.thickness;
    SPHSystem system(BoundingBoxd(Vecd::Constant(-2*R),Vecd::Constant(2*R)),o.spacing,o.threads);
    SolidBody core(system,makeShared<GeometricShapeBall>(Vecd::Zero(),rc,"Core"));
    core.defineMatterMaterial<DruckerPragerCapSolid>(o.core_rho,dpc::Model::fromCSV(o.dpc_table,o.dpc_r0,o.dpc_shear_dilation,o.dpc_shear_hardening));
    core.generateParticles<BaseParticles,CapsuleCoreGenerator>(rc,o.spacing);
    auto&cp=core.getBaseParticles();size_t cn=cp.TotalRealParticles();
    auto*cx=cp.getVariableDataByName<Vecd>("Position");auto*cv=cp.registerStateVariableData<Vecd>("Velocity");
    Field cX(cx,cx+cn);
    SnsniReference cr(cX,o.spacing,o.core_E/(2*(1+o.core_nu)),o.core_E*o.core_nu/((1+o.core_nu)*(1-2*o.core_nu)),o.core_rho,cells);
    SnsniPlasticSolid cb(cr,dpc::Model::fromCSV(o.dpc_table,o.dpc_r0,o.dpc_shear_dilation,o.dpc_shear_hardening));
    auto*cm=cp.getVariableDataByName<Real>("Mass");auto*cvol=cp.getVariableDataByName<Real>("VolumetricMeasure");
    for(size_t j=0;j<cn;++j){cm[j]=cr.mass[j];cvol[j]=cm[j]/o.core_rho;}

    PDBody shell(system,makeShared<DefaultShape>("NOSB_Shell"));
    shell.defineAdaptation<PDAdaptation>(1.5075,o.spacing/o.pd_spacing);
    shell.defineMatterMaterial<HughesWingetSolid>(o.shell_rho,o.shell_E,o.shell_nu);
    shell.generateParticles<BaseParticles,LayeredPDCapsuleGenerator>(rc,o.thickness,o.pd_spacing,o.shell_layers);
    InnerRelation relation(shell);
    SimpleDynamics<PDTimeStepInitialization> initialize(shell,makeShared<Gravity>(Vecd::Zero()));
    InteractionWithUpdate<NosbPDShapeMatrix> shape(relation);
    SimpleDynamics<NosbPDFirstStep> native_first(shell);
    InteractionWithUpdate<NosbPDSecondStep> native_second(relation);
    InteractionDynamics<NosbPDThirdStep> native_third(relation);
    InteractionDynamics<ProjectedPairStabilization> native_stabilizer(relation,o.pd_hourglass,.1);
    InteractionDynamics<PairNumericalDampingforPD> native_damping(relation);
    SimpleDynamics<RelaxedPDFourthStep> native_fourth(shell,o.damping_rate);
    auto&sp=shell.getBaseParticles();size_t sn=sp.TotalRealParticles(),columns_count=sn/o.shell_layers;
    auto*sx=sp.getVariableDataByName<Vecd>("Position");auto*sv=sp.getVariableDataByName<Vecd>("Velocity");
    auto*sm=sp.getVariableDataByName<Real>("Mass");auto*svol=sp.getVariableDataByName<Real>("VolumetricMeasure");
    Field sX(sx,sx+sn);StdVec<StdVec<size_t>>columns(columns_count);StdVec<Vecd>directors(columns_count);
    for(size_t c=0;c<columns_count;++c){for(size_t l=0;l<o.shell_layers;++l)columns[c].push_back(l*columns_count+c);directors[c]=sx[c].normalized();}
    ThicknessQuadrature rule(o.shell_layers);LayeredParticleGeometry geometry(sp,columns,directors,rule.xi);
    PDLayeredQuadraticReference sr(relation,geometry);
    InteractionWithUpdate<NosbPDQuadraticSecondStep> quadratic_second(relation,sr);
    InteractionDynamics<NosbPDQuadraticThirdStep> quadratic_third(relation,sr);
    PDEnergyCrackBand sb(relation,sr,o.shell_E,o.shell_nu,o.fracture_stress,o.fracture_energy,o.pd_hourglass,shell.getSPHAdaptation().getKernel()->CutOffRadius());
    system.initializeSystemCellLinkedLists();system.initializeSystemConfigurations();shape.exec();sr.build();sb.initialize();
    SolidForcePort core_port(core,SolidForcePort::Storage::ForcePrior),pd_port(shell,SolidForcePort::Storage::PDAccelerationPrior);
    StdVec<Vecd>points,outer_points;StdVec<Real>areas;
    for(const auto&face:cr.faces)if(face.boundary)for(const auto&q:face.q){points.push_back(q.x);areas.push_back(q.w);}
    for(const auto&n:sphereDirections(R,o.contact_quadrature_scale*o.spacing))outer_points.push_back(R*n);
    SnsniCoupledSurface cs(core_port,points,cr);
    RadialLayeredSolidSurface si(pd_port,points,geometry,rc,o.thickness,o.pd_spacing,o.interface_shell_support);
    LayeredSolidSurface so(pd_port,outer_points,geometry,1,o.pd_spacing,o.interface_shell_support);
    Real Eref=o.constraint_E>0?o.constraint_E:o.core_E,Ktie=o.tie_factor*Eref/o.spacing,Kcontact=o.contact_factor*Eref/o.spacing;
    BondedSolidInterface tie(si,cs,areas,Ktie,R);
    StdVec<Real>outer_areas(outer_points.size(),4*Pi*R*R/outer_points.size());
    SolidSurfacePlaneContact contact(so,outer_areas,Kcontact,R);
    size_t count=cn+sn;StdVec<Real>mass(count);Field X(count);
    for(size_t j=0;j<cn;++j){mass[j]=cm[j];X[j]=cX[j];}
    for(size_t j=0;j<sn;++j){mass[cn+j]=sm[j];X[cn+j]=sX[j];}
    SnsniCentralDifference time(mass);Field total(count,Vecd::Zero()),internal(count,Vecd::Zero());
    auto*sa=sp.getVariableDataByName<Vecd>("PDAcceleration");auto*sinternal=sp.getVariableDataByName<Vecd>("PDInternalForce");auto*shg=sp.getVariableDataByName<Vecd>("HourglassForce");
    auto*sF=sp.getVariableDataByName<Matd>("DeformationGradient");auto*sstress=sp.getVariableDataByName<Matd>("PDStress");
    Real core_mass=std::accumulate(cr.mass.begin(),cr.mass.end(),Real(0)),shell_mass=std::accumulate(sm,sm+sn,Real(0));
    Real exact_shell_mass=o.shell_rho*4*Pi*(R*R*R-rc*rc*rc)/3;
    if(std::abs(shell_mass/exact_shell_mass-1)>1e-10)throw std::runtime_error("Physical shell mass failed");
    Real radial_spacing=o.thickness;for(size_t k=1;k<rule.xi.size();++k)radial_spacing=std::min(radial_spacing,.5*o.thickness*(rule.xi[k]-rule.xi[k-1]));
    auto&mat=DynamicCast<HughesWingetSolid>(&shell,shell.getMatterMaterial());
    Real radial_dt=o.pd_cfl*std::min(radial_spacing,o.pd_spacing)/mat.ReferenceSoundSpeed();
    Real core_dt=o.cfl/std::sqrt(cb.linear_frequency_bound_squared);
    auto setPositions=[&](const Field&d,const Field&v){for(size_t j=0;j<cn;++j){cx[j]=cX[j]+d[j];cv[j]=v[j];}for(size_t j=0;j<sn;++j){sx[j]=sX[j]+d[cn+j];sv[j]=v[cn+j];}};
    auto coreEvaluate=[&](const Field&d,bool commit){Field u(d.begin(),d.begin()+cn);cb.evaluate(u,o.core_rho,o.cfl);if(commit)cb.commit();};
    auto shellEvaluate=[&](bool advance){tbb::parallel_for(size_t(0),sn,[&](size_t i){sb.evaluate(i,advance);});tbb::parallel_for(size_t(0),sn,[&](size_t i){sa[i].setZero();sb.forces(i);});};
    auto constraints=[&](Real compression){core_port.clear();pd_port.clear();tie.apply();contact.apply(-R+compression/2,R-compression/2);};
    auto forceAt=[&](const Field&d,const Field&v,Field&f,Real compression,bool advance){setPositions(d,v);coreEvaluate(d,true);shellEvaluate(advance);constraints(compression);f.resize(count);
        for(size_t j=0;j<cn;++j){internal[j]=cb.internal[j];f[j]=core_port.external[j]-internal[j];}
        for(size_t j=0;j<sn;++j){internal[cn+j]=-sinternal[j]-shg[j];f[cn+j]=sm[j]*pd_port.external[j]-internal[cn+j];}};
    std::ofstream parameters(out/"parameters.json");parameters<<std::setprecision(17)
      <<"{\"core_backend\":\"frozen_SNSNI_RKPM_DPC\",\"shell_backend\":\"NOSB_PD_layered_volume\",\"fluid_backend\":\"none\",\"shell_constitutive\":\"finite_volume_Biot_energy_damage\",\"fracture_representation\":\"diffuse_crack_band_not_sharp_crack\",\"interface\":\"bonded_vector_penalty_common_facet_trace\",\"plate_friction\":0,\"core_particles\":"<<cn<<",\"shell_particles\":"<<sn<<",\"material_sites\":"<<cr.cells.size()
      <<",\"core_diameter_m\":"<<2*rc<<",\"outer_diameter_m\":"<<2*R<<",\"shell_thickness_m\":"<<o.thickness<<",\"core_spacing_m\":"<<o.spacing<<",\"pd_spacing_m\":"<<o.pd_spacing
      <<",\"core_mass_kg\":"<<core_mass<<",\"shell_mass_kg\":"<<shell_mass<<",\"interface_K_Pa_per_m\":"<<Ktie<<",\"contact_K_Pa_per_m\":"<<Kcontact
      <<",\"shell_E_Pa\":"<<o.shell_E<<",\"shell_nu\":"<<o.shell_nu<<",\"ft_Pa\":"<<o.fracture_stress<<",\"Gc_J_m2\":"<<o.fracture_energy<<",\"crack_band_width_m\":"<<sb.fracture.width<<",\"kappa0\":"<<sb.fracture.k0<<",\"kappaf\":"<<sb.fracture.kf
      <<",\"core_linear_dt_s\":"<<core_dt<<",\"shell_initial_elastic_dt_s\":"<<sr.elastic_dt_bound<<",\"shell_radial_dt_s\":"<<radial_dt<<",\"shell_stabilizer_dt_s\":"<<sb.stabilizer_dt<<",\"interface_dt_s\":"<<tie.stable_dt<<",\"contact_dt_s\":"<<contact.stable_dt
      <<",\"core_affine_error\":"<<cs.affine_error<<",\"shell_inner_affine_error\":"<<si.affine_error<<",\"shell_outer_affine_error\":"<<so.affine_error<<",\"shell_polynomial_error\":"<<sr.max_reproduction_error<<",\"load_steps\":"<<o.load_steps<<",\"maximum_outer_strain\":"<<o.strain<<"}\n";parameters.close();

    if(o.preflight_only)
    {
        Field d(count),v(count);Matd affine=Matd::Identity();affine(0,0)=.998;affine(1,1)=1.001;affine(2,2)=.999;
        for(size_t j=0;j<count;++j){d[j]=(affine-Matd::Identity())*X[j];if(j<cn)d[j]+=Vecd(5e-8*rc/.0005,0,0);v[j]=Vecd(.2*std::sin(37*X[j][0]/R),-.1,.3*X[j][2]/R);}
        setPositions(d,v);coreEvaluate(d,false);shellEvaluate(false);
        Real affine_error=0,native_error=0;for(size_t l=0;l<cr.cells.size();++l){affine_error=std::max(affine_error,(cb.F[l]-affine).norm());auto a=cb.model.update(cb.F[l],cb.history[l]);native_error=std::max(native_error,(a.sigma-cb.pending[l].material.sigma).norm());}
        for(size_t j=0;j<sn;++j)affine_error=std::max(affine_error,(sF[j]-affine).norm());
        Real compression=.004*2*R;constraints(compression);Real power=0;
        for(size_t j=0;j<cn;++j)power+=core_port.external[j].dot(v[j]);for(size_t j=0;j<sn;++j)power+=sm[j]*pd_port.external[j].dot(v[cn+j]);
        Real eps=1e-12;Field plus=d,minus=d;for(size_t j=0;j<count;++j){plus[j]+=eps*v[j];minus[j]-=eps*v[j];}
        setPositions(plus,v);constraints(compression);Real ep=tie.energy+contact.energy;
        setPositions(minus,v);constraints(compression);Real em=tie.energy+contact.energy;
        Real derivative=(ep-em)/(2*eps),mapping_error=std::abs(derivative+power)/std::max({Real(1e-20),std::abs(derivative),std::abs(power)});
        cb.discard();auto stat=sb.statistics();
        bool pass=affine_error<1e-8&&native_error==0&&mapping_error<1e-5&&stat.initiated==0&&
          std::max({tie.worst_force,tie.worst_torque,tie.worst_power,contact.worst_force,contact.worst_torque,contact.worst_power})<1e-9;
        std::ofstream f(out/"preflight.json");f<<std::setprecision(17)<<"{\"status\":\""<<(pass?"PASSED":"FAILED")<<"\",\"affine_F_error\":"<<affine_error<<",\"native_DPC_stress_difference_Pa\":"<<native_error<<",\"constraint_energy_gradient_relative_error\":"<<mapping_error<<",\"interface_force_residual\":"<<tie.worst_force<<",\"interface_torque_residual\":"<<tie.worst_torque<<",\"interface_power_residual\":"<<tie.worst_power<<",\"contact_power_residual\":"<<contact.worst_power<<",\"advanced_core_histories\":0,\"advanced_shell_damage\":false,\"physical_histories\":0}\n";
        std::cout<<"PREFLIGHT "<<(pass?"PASSED":"FAILED")<<" affine="<<affine_error<<" mapping="<<mapping_error<<std::endl;return pass?0:2;
    }
    size_t stage=1;Real stage_time=0,elapsed=0;size_t checks=0;
    auto compression=[&](Real t){Real a=std::clamp(t/o.ramp_duration,Real(0),Real(1));return o.diameter*o.strain*(Real(stage-1)+.5*(1-std::cos(Pi*a)))/o.load_steps;};
    forceAt(time.displacement,time.velocity,total,0,false);for(size_t j=0;j<count;++j)time.acceleration[j]=total[j]/mass[j];
    std::ofstream progress(out/"progress.csv"),eq(out/"equilibria.csv"),events(out/"fracture_events.csv");
    const char*header="step,stage,elapsed_s,strain,compression_m,dt_s,upper_N,lower_N,force_residual,plate_imbalance,kinetic_ratio,core_min_J,shell_min_J,max_damage,initiated_sites,failed_sites,fracture_dissipation_J,shell_stored_J,interface_energy_J,interface_gap_m,contact_energy_J,max_penetration_m,platen_work_J,internal_work_J,damping_loss_J,budget_residual_J\n";
    progress<<std::setprecision(17)<<header;eq<<std::setprecision(17)<<header;events<<std::setprecision(17)<<"event,step,stage,strain,compression_m,mean_reaction_N,max_damage,initiated_sites,failed_sites\n";
    bool onset_saved=false,failure_saved=false;Real Wp=0,Wi=0,Wd=0,W0constraint=tie.energy+contact.energy,K0=time.kinetic();
    auto snapshot=[&](const std::string&tag,Real u){Field cpos(cr.cells.size()),cvelocity;StdVec<Real>p(cr.cells.size()),J(cr.cells.size()),alpha(cr.cells.size());Field cd(time.displacement.begin(),time.displacement.begin()+cn);
        for(size_t l=0;l<cr.cells.size();++l){cpos[l]=cr.cells[l].center+cr.value(l,cd);p[l]=cb.pending[l].material.p;J[l]=cb.F[l].determinant();alpha[l]=cb.history[l].alpha;}
        vtkPoints(out/"vtp"/("core_"+tag+".vtp"),cpos,{}, {},p,J,alpha,{},{});
        Field spts(sx,sx+sn),svel(time.velocity.begin()+cn,time.velocity.end()),normal(sb.crack_normal,sb.crack_normal+sn);StdVec<Real>sd(sb.damage,sb.damage+sn),spressure(sn),sJ(sn),loss(sb.dissipation_density,sb.dissipation_density+sn);
        for(size_t j=0;j<sn;++j){spressure[j]=-sstress[j].trace()/3;sJ[j]=sF[j].determinant();}
        vtkPoints(out/"vtp"/("shell_"+tag+".vtp"),spts,svel,sd,spressure,sJ,{},loss,normal);vtkPlates(out/"vtp"/("platens_"+tag+".vtp"),R,u);};
    if(recovery::checkpoint.resume.empty())snapshot("0",0);
    auto transfer=[&](recovery::Archive&a){
        a.fields(stage,stage_time,elapsed,checks,onset_saved,failure_saved,Wp,Wi,Wd,W0constraint,K0,total,internal);
        recovery::integrator(a,time);recovery::solid(a,cb);recovery::contact(a,contact);
        a.fields(tie.energy,tie.max_gap,tie.force_residual,tie.torque_residual,tie.power_residual,tie.worst_force,tie.worst_torque,tie.worst_power);
        a.value(sb.history);a.span(sb.damage,sn);a.span(sb.dissipation_density,sn);a.span(sb.stored_density,sn);
        a.span(sb.hg_density,sn);a.span(sb.kappa_field,sn);a.span(sb.tangent_bound,sn);a.span(sb.material_failed,sn);a.span(sb.crack_normal,sn);
        a.span(sF,sn);a.span(sstress,sn);a.span(sa,sn);a.span(sinternal,sn);a.span(shg,sn);
        a.span(sp.getVariableDataByName<Matd>("PDPK1"),sn);
        a.span(cx,cn);a.span(cv,cn);a.span(sx,sn);a.span(sv,sn);
        if(time.displacement.size()!=count||cb.history.size()!=cr.cells.size()||sb.history.size()!=sn)throw std::runtime_error("Checkpoint geometry size mismatch");
    };
    auto save=[&](){progress.flush();eq.flush();events.flush();recovery::checkpoint.save(out,{"progress.csv","equilibria.csv","fracture_events.csv"},time.steps,transfer);};
    if(!recovery::checkpoint.resume.empty()){
        progress.close();eq.close();events.close();
        recovery::checkpoint.load(out,{"progress.csv","equilibria.csv","fracture_events.csv"},transfer);
        progress.open(out/"progress.csv",std::ios::app);eq.open(out/"equilibria.csv",std::ios::app);events.open(out/"fracture_events.csv",std::ios::app);
        progress<<std::setprecision(17);eq<<std::setprecision(17);events<<std::setprecision(17);
    }
    save();
    while(stage<=o.load_steps)
    {
        size_t prior_stage=stage;Real shell_dt=sb.finiteStrainTimeStep(sr.elastic_dt_bound);
        Real bound=1/std::sqrt(1/(core_dt*core_dt)+1/(shell_dt*shell_dt)+1/(sb.stabilizer_dt*sb.stabilizer_dt)+1/(tie.stable_dt*tie.stable_dt)+1/(contact.stable_dt*contact.stable_dt));
        Real dt=std::min({o.max_dt,bound,cb.acoustic_dt,radial_dt});if(stage_time<o.ramp_duration)dt=std::min(dt,o.ramp_duration-stage_time);
        if(!std::isfinite(dt)||dt<1e-15||time.steps>=o.max_steps)throw std::runtime_error("Coupled time/count bound exhausted; no retry or force correction");
        Field old_d=time.displacement,old_internal=internal;Real old_force=.5*(contact.upper+contact.lower),old_u=compression(stage_time),oldK=time.kinetic();
        time.step(dt,o.damping_rate,[&](const Field&d,const Field&v,Field&f){forceAt(d,v,f,compression(stage_time+dt),true);total=f;});
        Real u=compression(stage_time+dt),mean=.5*(contact.upper+contact.lower),newK=time.kinetic(),factor2=std::exp(-o.damping_rate*dt);
        Wp+=.5*(old_force+mean)*(u-old_u);for(size_t j=0;j<count;++j)Wi+=.5*(old_internal[j]+internal[j]).dot(time.displacement[j]-old_d[j]);
        Wd+=oldK*(1-factor2)+newK*(1/factor2-1);stage_time+=dt;elapsed+=dt;
        auto stat=sb.statistics();
        if((!onset_saved&&stat.initiated)||(!failure_saved&&stat.failed))
        {std::string event=!onset_saved?"first_damage":"first_complete_tensile_failure";events<<event<<','<<time.steps<<','<<stage<<','<<u/o.diameter<<','<<u<<','<<mean<<','<<stat.max_damage<<','<<stat.initiated<<','<<stat.failed<<'\n';events.flush();snapshot(event,u);if(!onset_saved)onset_saved=true;else failure_saved=true;}
        if(time.steps%200==0)
        {
            Real residual=0,minc=1,mins=1;for(const auto&f:total)residual+=f.squaredNorm();residual=std::sqrt(residual)/std::max(mean,Real(1e-20));
            for(const auto&F:cb.F)minc=std::min(minc,F.determinant());for(size_t j=0;j<sn;++j)mins=std::min(mins,sF[j].determinant());
            Real imbalance=std::abs(contact.upper-contact.lower)/std::max(mean,Real(1e-20)),kr=newK/std::max(mean*u,Real(1e-30));
            Real budget=Wp-Wi-(tie.energy+contact.energy-W0constraint)-(newK-K0)-Wd;
            auto row=[&](std::ostream&f){f<<time.steps<<','<<stage<<','<<elapsed<<','<<u/o.diameter<<','<<u<<','<<dt<<','<<contact.upper<<','<<contact.lower<<','<<residual<<','<<imbalance<<','<<kr<<','<<minc<<','<<mins<<','<<stat.max_damage<<','<<stat.initiated<<','<<stat.failed<<','<<stat.dissipated<<','<<stat.stored<<','<<tie.energy<<','<<tie.max_gap<<','<<contact.energy<<','<<contact.max_penetration<<','<<Wp<<','<<Wi<<','<<Wd<<','<<budget<<'\n';f.flush();};
            row(progress);if(time.steps%1000==0)std::cout<<"PROGRESS step="<<time.steps<<" stage="<<stage<<" strain="<<u/o.diameter<<" force="<<mean<<" residual="<<residual<<" KE_ratio="<<kr<<" damage="<<stat.max_damage<<" dt="<<dt<<std::endl;
            bool good=stage_time>=o.ramp_duration&&mean>1e-10&&residual<=o.force_tolerance&&imbalance<=o.balance_tolerance&&kr<=o.kinetic_tolerance;checks=good?checks+1:0;
            if(checks>=3){row(eq);snapshot(std::to_string(stage),u);std::cout<<"ACCEPTED stage="<<stage<<" force="<<mean<<" damage="<<stat.max_damage<<std::endl;++stage;stage_time=0;checks=0;}
        }
        if(stage_time>o.max_stage_duration)throw std::runtime_error("Coupled stage did not equilibrate; no completed curve claimed");
        if(time.steps%200==0 && recovery::checkpoint.due(time.steps,stage!=prior_stage))save();
        if(recovery::checkpoint.stop_step && time.steps>=recovery::checkpoint.stop_step){save();return 86;}
    }
    std::ofstream summary(out/"summary.json");summary<<"{\"status\":\"COMPLETE_FROZEN_CORE_NOSB_SNSNI_CAPSULE\",\"accepted_points\":"<<o.load_steps<<",\"steps\":"<<time.steps<<",\"core_material_changed\":false,\"parameters_fitted\":false,\"observed_damage_onset\":"<<(onset_saved?"true":"false")<<",\"observed_complete_tensile_failure\":"<<(failure_saved?"true":"false")<<"}\n";
    return 0;
}

int main(int argc,char**argv)
{
    try{std::string cells;std::vector<char*>args{argv[0]};for(int i=1;i<argc;++i){std::string a=argv[i];if(a.rfind("--snsni-cells=",0)==0)cells=a.substr(14);else if(recovery::checkpoint.take(a)){}else args.push_back(argv[i]);}
        if(cells.empty())throw std::invalid_argument("--snsni-cells is required");HybridOptions options;options.parseHybrid(int(args.size()),args.data());
        tbb::task_arena arena(int(options.threads));return arena.execute([&]{return runSnsniCapsule(options,cells);});
    }catch(const std::exception&e){std::cerr<<"COUPLED_ERROR "<<e.what()<<std::endl;return 1;}
}
