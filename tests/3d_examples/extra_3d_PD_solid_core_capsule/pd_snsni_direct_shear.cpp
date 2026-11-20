/** S-NSNI/RKPM DPC core + layered NOSB-PD shell with optional transport exports.
 * One common explicit displacement/velocity state, physical lumped masses,
 * common-trace vector penalty and frictionless reconstructed outer contact.
 * No WCSPH fluid, material retuning, particle erosion or force rescaling.
 */
#include "pd_snsni_capsule_fixture.h"
#include "solid_snsni_plastic.h"
#include "pd_radial_layered_surface.h"
#include "recovery_checkpoint.h"
#include "diffuse_crack_release.h"
#include "rigid_grain_shear.h"
#if defined(__GLIBC__)
#include <malloc.h>
#endif
using namespace SPH::solid_dynamics;
using namespace SPH::peridynamics;
using Field=StdVec<Vecd>;
std::string snapshot_barrier_input;

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

int runDirectShear(const HybridOptions&o,const std::string&cells,const std::string&sand,const ShearBoxOptions&box,const PeriExecution&execution)
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
    SnsniPlasticSolid cb(cr,dpc::Model::fromCSV(o.dpc_table,o.dpc_r0,o.dpc_shear_dilation,o.dpc_shear_hardening),
                         execution);
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
    if(!snapshot_barrier_input.empty())
    {
        // Kinematics-only reconstruction from an actual native VTP snapshot.
        // No core or shell history update, force evaluation or time step.
        std::ifstream stream(snapshot_barrier_input);if(!stream)throw std::runtime_error("Missing shell snapshot");
        std::string xml((std::istreambuf_iterator<char>(stream)),{});
        auto values=[&](const std::string&name,bool points=false){
            auto start=xml.find(points?"<Points>":("Name=\""+name+"\""));
            if(start==std::string::npos)throw std::runtime_error("Missing VTP field: "+name);
            if(points)start=xml.find("<DataArray",start);
            auto begin=xml.find('>',start),end=xml.find("</DataArray>",begin);
            if(begin==std::string::npos||end==std::string::npos)throw std::runtime_error("Incomplete VTP field");
            std::istringstream in(xml.substr(begin+1,end-begin-1));std::vector<Real>result;Real x;
            while(in>>x)result.push_back(x);return result;};
        auto pos=values("Position",true),d=values("Damage"),n=values("CrackNormalReference"),j=values("J");
        if(pos.size()!=3*sn||d.size()!=sn||n.size()!=3*sn||j.size()!=sn)throw std::runtime_error("VTP shell count mismatch");
        auto*f=sp.getVariableDataByName<Matd>("DeformationGradient");
        for(size_t i=0;i<sn;++i){sx[i]=Vecd(pos[3*i],pos[3*i+1],pos[3*i+2]);sb.damage[i]=d[i];sb.crack_normal[i]=Vecd(n[3*i],n[3*i+1],n[3*i+2]);}
        Real error=0;
        for(size_t i=0;i<sn;++i){f[i].setZero();for(const auto&e:sr.rows[i])f[i]+=(sx[e.j]-sx[i])*e.gradient.transpose();error=std::max(error,std::abs(f[i].determinant()-j[i]));}
        if(error>1e-8)throw std::runtime_error("Rebuilt snapshot F does not reproduce stored J");
        writeLayeredReleaseBarrier((out/"release_barrier.csv").string(),geometry,sb,f,directors,rule.weight);
        std::ofstream meta(out/"rebuild.json");meta<<std::setprecision(17)<<"{\"physical_histories_advanced\":0,\"stored_J_max_error\":"<<error<<",\"source\":"<<std::quoted(snapshot_barrier_input)<<"}\n";
        return 0;
    }

    SolidForcePort core_port(core,SolidForcePort::Storage::ForcePrior),pd_port(shell,SolidForcePort::Storage::PDAccelerationPrior);
    StdVec<Vecd>points,outer_points;StdVec<Real>areas;
    for(const auto&face:cr.faces)if(face.boundary)for(const auto&q:face.q){points.push_back(q.x);areas.push_back(q.w);}
    for(const auto&n:sphereDirections(R,o.contact_quadrature_scale*o.spacing))outer_points.push_back(R*n);
    SnsniCoupledSurface cs(core_port,points,cr);
    RadialLayeredSolidSurface si(pd_port,points,geometry,rc,o.thickness,o.pd_spacing,o.interface_shell_support);
    LayeredSolidSurface so(pd_port,outer_points,geometry,1,o.pd_spacing,o.interface_shell_support);
    Real Eref=o.constraint_E>0?o.constraint_E:o.core_E,Ktie=o.tie_factor*Eref/o.spacing,Kcontact=o.contact_factor*Eref/o.spacing;
    BondedSolidInterface tie(si,cs,areas,Ktie,R,execution.makeBondedInterface());
    StdVec<Real>outer_areas(outer_points.size(),4*Pi*R*R/outer_points.size());
    RigidGrainShearBox world(so,outer_areas,Kcontact,sand,box,&cs,areas);
    const size_t released_preprocessing_bytes=cr.releasePreprocessingCache();
#if defined(__GLIBC__)
    // Return the one-time construction buffers to the OS, outside time stepping.
    malloc_trim(0);
#endif
    world.cacheCorePositions(tie.x_b);
    std::cout<<"OPTIMIZED released_preprocessing_bytes="<<released_preprocessing_bytes<<std::endl;
    const size_t ncaps=cn+sn,ng=world.grains.size(),count=ncaps+2*ng;
    StdVec<Real>mass(count);Field X(count,Vecd::Zero());
    for(size_t j=0;j<cn;++j){mass[j]=cm[j];X[j]=cX[j];}
    for(size_t j=0;j<sn;++j){mass[cn+j]=sm[j];X[cn+j]=sX[j];}
    for(size_t j=0;j<ng;++j){mass[ncaps+2*j]=world.grains[j].mass;mass[ncaps+2*j+1]=world.grains[j].inertia;X[ncaps+2*j]=world.grains[j].reference;}
    SnsniCentralDifference time(mass);Field total(count,Vecd::Zero()),internal(count,Vecd::Zero());
    auto*sa=sp.getVariableDataByName<Vecd>("PDAcceleration");auto*sinternal=sp.getVariableDataByName<Vecd>("PDInternalForce");auto*shg=sp.getVariableDataByName<Vecd>("HourglassForce");
    auto*sF=sp.getVariableDataByName<Matd>("DeformationGradient");auto*sstress=sp.getVariableDataByName<Matd>("PDStress");
    auto&mat=DynamicCast<HughesWingetSolid>(&shell,shell.getMatterMaterial());
    Real radial_spacing=o.thickness;for(size_t k=1;k<rule.xi.size();++k)radial_spacing=std::min(radial_spacing,.5*o.thickness*(rule.xi[k]-rule.xi[k-1]));
    Real radial_dt=o.pd_cfl*std::min(radial_spacing,o.pd_spacing)/mat.ReferenceSoundSpeed();
    Real core_dt=o.cfl/std::sqrt(cb.linear_frequency_bound_squared);
    auto setPositions=[&](const Field&d,const Field&v){
        for(size_t j=0;j<cn;++j){cx[j]=cX[j]+d[j];cv[j]=v[j];}
        for(size_t j=0;j<sn;++j){sx[j]=sX[j]+d[cn+j];sv[j]=v[cn+j];}
        for(size_t j=0;j<ng;++j){auto&g=world.grains[j];g.x=g.reference+d[ncaps+2*j];g.v=v[ncaps+2*j];g.omega=v[ncaps+2*j+1];}
    };
    auto constraints=[&](Field&f,Real shear,Real speed,Real load,Real dt,bool commit){
        core_port.clear();pd_port.clear();tie.apply();Field bonded(sn),bonded_core(cn);for(size_t j=0;j<sn;++j)bonded[j]=sm[j]*pd_port.external[j];for(size_t j=0;j<cn;++j)bonded_core[j]=core_port.external[j];
        core_port.clear();pd_port.clear();world.evaluate(shear,speed,load,dt,commit);f.assign(count,Vecd::Zero());
        for(size_t j=0;j<cn;++j)f[j]=bonded_core[j]+core_port.external[j];
        for(size_t j=0;j<sn;++j)f[cn+j]=bonded[j]+sm[j]*pd_port.external[j];
        for(size_t j=0;j<ng;++j){f[ncaps+2*j]=world.grains[j].force;f[ncaps+2*j+1]=world.grains[j].torque;}
    };
    auto forceAt=[&](const Field&d,const Field&v,Field&f,Real shear,Real speed,Real load,Real dt,bool advance){
        setPositions(d,v);Field u(d.begin(),d.begin()+cn);cb.evaluate(u,o.core_rho,o.cfl);cb.commit();
        tbb::parallel_for(size_t(0),sn,[&](size_t i){sb.evaluate(i,advance);});
        tbb::parallel_for(size_t(0),sn,[&](size_t i){sa[i].setZero();sb.forces(i);});
        constraints(f,shear,speed,load,dt,advance);internal.assign(count,Vecd::Zero());
        for(size_t j=0;j<cn;++j){internal[j]=cb.internal[j];f[j]-=internal[j];}
        for(size_t j=0;j<sn;++j){internal[cn+j]=-sinternal[j]-shg[j];f[cn+j]-=internal[cn+j];}
    };
    Real capsule_mass=std::accumulate(cm,cm+cn,Real(0))+std::accumulate(sm,sm+sn,Real(0)),grain_mass=0;
    for(auto&g:world.grains)grain_mass+=g.mass;
    std::ofstream parameters(out/"parameters.json");parameters<<std::setprecision(17)
      <<"{\"core_backend\":\"S-NSNI_RKPM_DPC\",\"shell_backend\":\"NOSB_PD_layered_volume\",\"sand_backend\":\"rotating_rigid_spheres_Hertz_linear_history_Coulomb\",\"core_particles\":"<<cn<<",\"shell_particles\":"<<sn<<",\"material_sites\":"<<cr.cells.size()<<",\"sand_grains\":"<<ng
      <<",\"execution_backend\":"<<std::quoted(execution.backendName())
      <<",\"core_derivatives_backend\":"<<std::quoted(cb.hasDeviceDerivatives()?"cuda_fp64":"cpu_fp64")
      <<",\"core_internal_force_backend\":"<<std::quoted(cb.hasDeviceDerivatives()?"cuda_fp64":"cpu_fp64")
      <<",\"core_material_backend\":"<<std::quoted(cb.hasDeviceMaterial()?"cuda_fp64":"cpu_fp64")<<",\"trial_stress_gradient_backend\":"<<std::quoted(cb.hasDeviceMaterial()?"cuda_fp64":"cpu_fp64")
      <<",\"bonded_interface_backend\":"<<std::quoted(execution.bondedInterfaceOnDevice()?"cuda_fp64":"cpu_fp64")
      <<",\"capsule_count\":1,\"capsule_mass_kg\":"<<capsule_mass<<",\"sand_mass_kg\":"<<grain_mass<<",\"capsule_mass_fraction\":"<<capsule_mass/(grain_mass+capsule_mass)<<",\"capsule_sand_mu\":0,\"initial_height_m\":"<<box.initial_height_m
      <<",\"initial_area_m2\":"<<box.length_x_m*box.length_y_m<<",\"normal_pressure_Pa\":"<<box.normal_pressure_Pa<<",\"normal_load_control\":\"massless_exact_piston_constant_nominal_force\",\"load_steps\":"<<box.normal_increments+box.shear_increments
      <<",\"core_material_changed\":false,\"shell_material_changed\":false,\"pore_pressure_solved\":false,\"chemical_strength_feedback\":false,\"representative_volume_claimed\":false,\"shear_rate_validated\":false}\n";parameters.close();

    if(o.preflight_only){
        Field d(count,Vecd::Zero()),v(count,Vecd::Zero()),f;
        // Exercise an ACTIVE sphere/capsule contact, rather than accepting a
        // vacuous zero-force transpose check on the barely touching packing.
        size_t nearest=0;Real smallest=1e10;
        for(size_t i=0;i<ng;++i){const auto&g=world.grains[i];Real gap=g.reference.norm()-g.radius-R;if(gap<smallest){smallest=gap;nearest=i;}}
        d[ncaps+2*nearest]=-(smallest+1e-6)*world.grains[nearest].reference.normalized();
        for(size_t i=0;i<count;++i)v[i]=Vecd(std::sin(Real(i+1)),.3*std::cos(Real(2*i+1)),.1);
        auto potential=[&](const Field&u){setPositions(u,v);constraints(f,0,0,box.normal_pressure_Pa*box.length_x_m*box.length_y_m,0,false);
            return tie.energy+world.elastic_energy+world.normal_load*(world.top-box.initial_height_m/2);};
        potential(d);Real power=0;for(size_t i=0;i<count;++i)power+=f[i].dot(v[i]);
        const Real eps=1e-13;Field plus=d,minus=d;for(size_t i=0;i<count;++i){plus[i]+=eps*v[i];minus[i]-=eps*v[i];}
        Real ep=potential(plus),em=potential(minus),derivative=(ep-em)/(2*eps);
        Real error=std::abs(derivative+power)/std::max({Real(1e-20),std::abs(derivative),std::abs(power)});
        Real fn=world.hertz(1e-7,1e-4);Real expected=Real(4)/3*(box.grain_E_Pa/(2*(1-box.grain_nu*box.grain_nu)))*1e-2*std::pow(1e-7,1.5);
        ObjectiveLinearFriction history;history.force(Vecd::UnitZ(),Vecd::Zero(),0,10,1000,.5,0,true);
        Vecd ft=history.force(Vecd::UnitZ(),Vecd(1,0,0),0,10,1000,.5,.01,true);
        Real friction_error=std::abs(ft.norm()-5),before=history.stored;
        Vecd rotated=history.force(Vecd::UnitX(),Vecd::Zero(),0,10,1000,.5,0,true);
        Real rotation_error=std::abs(history.stored-before);
        bool pass=world.capsule_contacts>0&&error<1e-5&&std::abs(fn/expected-1)<1e-12&&friction_error<1e-12&&rotation_error<1e-12&&history.dissipated>=0&&world.normal_control_residual<1e-10&&std::max(world.force_mapping_residual,world.torque_mapping_residual)<1e-10;
        std::ofstream gate(out/"preflight.json");gate<<std::setprecision(17)<<"{\"status\":\""<<(pass?"PASSED":"FAILED")<<"\",\"active_capsule_contact_quadrature_points\":"<<world.capsule_contacts<<",\"constraint_potential_gradient_relative_error\":"<<error<<",\"Hertz_relative_error\":"<<std::abs(fn/expected-1)<<",\"Coulomb_force_error_N\":"<<friction_error<<",\"contact_frame_rotation_energy_error_J\":"<<rotation_error<<",\"normal_control_residual\":"<<world.normal_control_residual<<",\"force_mapping_residual\":"<<world.force_mapping_residual<<",\"torque_mapping_residual\":"<<world.torque_mapping_residual<<",\"core_histories_advanced\":0,\"shell_damage_advanced\":false,\"shear_curve_generated\":false}\n";
        return pass?0:2;
    }
    size_t stage=1,checks=0;Real stage_time=0,elapsed=0;bool damage_saved=false,failure_saved=false,access_saved=false;
    struct Load{Real shear,speed,normal;};
    auto loading=[&](Real t){Real duration=stage<=box.normal_increments?box.normal_ramp_s:box.shear_ramp_s;
        Real a=std::clamp(t/duration,Real(0),Real(1)),r=.5*(1-std::cos(Pi*a)),dr=t<duration?.5*Pi/duration*std::sin(Pi*a):0;
        if(stage<=box.normal_increments)return Load{0,0,box.normal_pressure_Pa*box.length_x_m*box.length_y_m*(Real(stage-1)+r)/box.normal_increments};
        Real dx=box.shear_end_m/box.shear_increments;return Load{dx*(Real(stage-box.normal_increments-1)+r),dx*dr,box.normal_pressure_Pa*box.length_x_m*box.length_y_m};};
    forceAt(time.displacement,time.velocity,total,0,0,0,0,false);for(size_t j=0;j<count;++j)time.acceleration[j]=total[j]/mass[j];
    Real Wp=0,Wi=0,Wd=0,W0=tie.energy+world.elastic_energy+world.tangent_energy,K0=time.kinetic();
    std::ofstream progress(out/"progress.csv"),eq(out/"equilibria.csv"),events(out/"fracture_events.csv");
    const char*header="step,stage,elapsed_s,shear_m,normal_load_N,shear_force_N,shear_stress_Pa,normal_stress_Pa,stress_ratio,top_height_m,dt_s,force_residual,imbalance,kinetic_ratio,core_min_J,shell_min_J,max_damage,initiated_sites,failed_sites,fracture_dissipation_J,interface_energy_J,contact_energy_J,friction_dissipation_J,interface_gap_m,max_penetration_m,max_overlap_ratio,grain_contacts,capsule_contacts,platen_work_J,internal_work_J,damping_loss_J,budget_residual_J\n";
    progress<<std::setprecision(17)<<header;eq<<std::setprecision(17)<<header;
    events<<std::setprecision(17)<<"event,step,stage,shear_m,normal_load_N,shear_force_N,max_damage,initiated_sites,failed_sites\n";
    auto snapshot=[&](const std::string&tag){
        Field cpos(cr.cells.size()),cd(time.displacement.begin(),time.displacement.begin()+cn);StdVec<Real>p(cr.cells.size()),J(cr.cells.size()),alpha(cr.cells.size());
        for(size_t l=0;l<cr.cells.size();++l){cpos[l]=cr.cells[l].center+cr.value(l,cd);p[l]=cb.pending[l].material.p;J[l]=cb.F[l].determinant();alpha[l]=cb.history[l].alpha;}
        vtkPoints(out/"vtp"/("core_"+tag+".vtp"),cpos,{},{},p,J,alpha,{},{});
        Field spts(sx,sx+sn),svel(time.velocity.begin()+cn,time.velocity.begin()+ncaps),normal(sb.crack_normal,sb.crack_normal+sn);StdVec<Real>sd(sb.damage,sb.damage+sn),spressure(sn),sJ(sn),loss(sb.dissipation_density,sb.dissipation_density+sn);
        for(size_t j=0;j<sn;++j){spressure[j]=-sstress[j].trace()/3;sJ[j]=sF[j].determinant();}
        vtkPoints(out/"vtp"/("shell_"+tag+".vtp"),spts,svel,sd,spressure,sJ,{},loss,normal);
        std::ofstream sandvtk(out/"vtp"/("sand_"+tag+".vtp"));sandvtk<<std::setprecision(17)<<"<VTKFile type=\"PolyData\" version=\"0.1\"><PolyData><Piece NumberOfPoints=\""<<ng<<"\" NumberOfVerts=\""<<ng<<"\"><Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
        for(auto&g:world.grains)sandvtk<<g.x[0]<<' '<<g.x[1]<<' '<<g.x[2]<<' ';sandvtk<<"</DataArray></Points><Verts><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">";for(size_t j=0;j<ng;++j)sandvtk<<j<<' ';sandvtk<<"</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">";for(size_t j=0;j<ng;++j)sandvtk<<j+1<<' ';sandvtk<<"</DataArray></Verts><PointData><DataArray type=\"Float64\" Name=\"Radius_m\" format=\"ascii\">";for(auto&g:world.grains)sandvtk<<g.radius<<' ';sandvtk<<"</DataArray></PointData></Piece></PolyData></VTKFile>\n";
        writeLayeredReleaseBarrier((out/("release_barrier_"+tag+".csv")).string(),geometry,sb,sF,directors,rule.weight);
        std::ofstream cellsout(out/("transport_cells_"+tag+".csv"));cellsout<<std::setprecision(17)<<"cell,J,relative_density,F00,F01,F02,F10,F11,F12,F20,F21,F22\n";
        for(size_t l=0;l<cr.cells.size();++l){cellsout<<l<<','<<cb.F[l].determinant()<<','<<cb.pending[l].material.relative_density;for(int a=0;a<3;++a)for(int b=0;b<3;++b)cellsout<<','<<cb.F[l](a,b);cellsout<<'\n';}
        std::ofstream state(out/("snapshot_"+tag+".json"));auto load=loading(stage_time);state<<std::setprecision(17)<<"{\"stage\":"<<stage<<",\"shear_m\":"<<load.shear<<",\"normal_load_N\":"<<load.normal<<",\"top_height_m\":"<<world.top<<",\"sand_geometry_is_explicit\":true,\"geometric_solid_extrusion_claimed\":false}\n";
    };
    auto transfer=[&](recovery::Archive&a){
        a.fields(stage,stage_time,elapsed,checks,damage_saved,failure_saved,access_saved,Wp,Wi,Wd,W0,K0,total,internal);
        recovery::integrator(a,time);recovery::solid(a,cb);
        a.fields(tie.energy,tie.max_gap,tie.force_residual,tie.torque_residual,tie.power_residual,tie.worst_force,tie.worst_torque,tie.worst_power);
        a.value(sb.history);a.span(sb.damage,sn);a.span(sb.dissipation_density,sn);a.span(sb.stored_density,sn);a.span(sb.hg_density,sn);a.span(sb.kappa_field,sn);a.span(sb.tangent_bound,sn);a.span(sb.material_failed,sn);a.span(sb.crack_normal,sn);
        a.span(sF,sn);a.span(sstress,sn);a.span(sa,sn);a.span(sinternal,sn);a.span(shg,sn);a.span(sp.getVariableDataByName<Matd>("PDPK1"),sn);a.span(cx,cn);a.span(cv,cn);a.span(sx,sn);a.span(sv,sn);
        a.fields(world.top,world.normal_load,world.elastic_energy,world.tangent_energy,world.friction_dissipation,world.max_sphere_penetration,world.max_capsule_penetration,world.max_overlap_ratio,world.stable_dt,world.lower_reaction,world.upper_reaction,world.boundary_force,world.boundary_torque);
        for(auto&g:world.grains)a.fields(g.x,g.v,g.omega,g.force,g.torque);
        for(auto&h:world.pairs)a.fields(h.slip,h.normal,h.active,h.stored,h.dissipated,h.released);
        for(auto&h:world.walls)a.fields(h.slip,h.normal,h.active,h.stored,h.dissipated,h.released);
        if(a.reading)world.invalidateNeighborCache();
    };
    auto save=[&](){progress.flush();eq.flush();events.flush();recovery::checkpoint.save(out,{"progress.csv","equilibria.csv","fracture_events.csv"},time.steps,transfer);};
    if(!recovery::checkpoint.resume.empty()){
        progress.close();eq.close();events.close();recovery::checkpoint.load(out,{"progress.csv","equilibria.csv","fracture_events.csv"},transfer);
        progress.open(out/"progress.csv",std::ios::app);eq.open(out/"equilibria.csv",std::ios::app);events.open(out/"fracture_events.csv",std::ios::app);progress<<std::setprecision(17);eq<<std::setprecision(17);events<<std::setprecision(17);
    }else snapshot("initial");
    save();auto lastsave=std::chrono::steady_clock::now();
    while(stage<=box.normal_increments+box.shear_increments){
        size_t prior_stage=stage;Real ramp=stage<=box.normal_increments?box.normal_ramp_s:box.shear_ramp_s;
        Real shell_dt=sb.finiteStrainTimeStep(sr.elastic_dt_bound);
        Real bound=1/std::sqrt(1/(core_dt*core_dt)+1/(shell_dt*shell_dt)+1/(sb.stabilizer_dt*sb.stabilizer_dt)+1/(tie.stable_dt*tie.stable_dt)+1/(world.stable_dt*world.stable_dt));
        Real dt=std::min({o.max_dt,bound,cb.acoustic_dt,radial_dt});if(stage_time<ramp)dt=std::min(dt,ramp-stage_time);
        if(!std::isfinite(dt)||dt<1e-15||time.steps>=o.max_steps)throw std::runtime_error("Direct-shear coupled step/count bound exhausted");
        Field old_d=time.displacement,old_internal=internal;Real old_force=-world.lower_reaction[0],old_height=world.top,oldK=time.kinetic();Load old=loading(stage_time);
        Load next=loading(stage_time+dt);
        time.step(dt,o.damping_rate,[&](const Field&d,const Field&v,Field&f){forceAt(d,v,f,next.shear,next.speed,next.normal,dt,true);total=f;});
        Real newK=time.kinetic(),factor2=std::exp(-o.damping_rate*dt),shear_force=-world.lower_reaction[0];
        Wp+=.5*(old_force+shear_force)*(next.shear-old.shear)-.5*(old.normal+next.normal)*(world.top-old_height);
        for(size_t j=0;j<ncaps;++j)Wi+=.5*(old_internal[j]+internal[j]).dot(time.displacement[j]-old_d[j]);
        Wd+=oldK*(1-factor2)+newK*(1/factor2-1);stage_time+=dt;elapsed+=dt;auto stat=sb.statistics();
        auto event=[&](const std::string&name){events<<name<<','<<time.steps<<','<<stage<<','<<next.shear<<','<<next.normal<<','<<shear_force<<','<<stat.max_damage<<','<<stat.initiated<<','<<stat.failed<<'\n';events.flush();snapshot(name);};
        if(!damage_saved&&stat.initiated){event("first_damage");damage_saved=true;}
        if(!failure_saved&&stat.failed){event("first_complete_tensile_failure");failure_saved=true;}
        if(!access_saved&&stat.initiated&&time.steps%200==0){bool connected=false;for(size_t c=0;c<geometry.columns.size()&&!connected;++c){StdVec<Matd>F;StdVec<Real>d;StdVec<Vecd>n;bool all=true;
            for(size_t j:geometry.columns[c]){all=all&&sb.damage[j]>0;F.push_back(sF[j]);d.push_back(sb.damage[j]);n.push_back(sb.crack_normal[j]);}
            if(all)connected=diffuseLayerBarrier(F,d,n,rule.weight,directors[c],geometry.thickness[c],sb.fracture.width).open_fraction>0;}
            if(connected){event("first_hydraulic_access");access_saved=true;}}
        if(time.steps%box.output_interval_steps==0){
            Real residual=0,minc=1,mins=1;for(size_t j=0;j<count;++j)residual+=total[j].squaredNorm()/(j>=ncaps&&(j-ncaps)%2?std::pow(world.grains[(j-ncaps)/2].radius,2):1);
            residual=std::sqrt(residual)/std::max(next.normal,Real(1e-20));
            for(const auto&F:cb.F)minc=std::min(minc,F.determinant());for(size_t j=0;j<sn;++j)mins=std::min(mins,sF[j].determinant());
            Real imbalance=(world.lower_reaction+world.upper_reaction).norm()/std::max(next.normal,Real(1e-20));
            Real kr=newK/std::max(next.normal*Real(.241e-3),Real(1e-30)),budget=Wp-Wi-(tie.energy+world.elastic_energy+world.tangent_energy-W0)-(newK-K0)-Wd-world.friction_dissipation;
            auto row=[&](std::ostream&f){Real A=box.length_x_m*box.length_y_m;f<<time.steps<<','<<stage<<','<<elapsed<<','<<next.shear<<','<<next.normal<<','<<shear_force<<','<<shear_force/A<<','<<next.normal/A<<','<<shear_force/std::max(next.normal,Real(1e-30))<<','<<world.top<<','<<dt<<','<<residual<<','<<imbalance<<','<<kr<<','<<minc<<','<<mins<<','<<stat.max_damage<<','<<stat.initiated<<','<<stat.failed<<','<<stat.dissipated<<','<<tie.energy<<','<<world.elastic_energy+world.tangent_energy<<','<<world.friction_dissipation<<','<<tie.max_gap<<','<<std::max(world.max_sphere_penetration,world.max_capsule_penetration)<<','<<world.max_overlap_ratio<<','<<world.grain_contacts<<','<<world.capsule_contacts<<','<<Wp<<','<<Wi<<','<<Wd<<','<<budget<<'\n';f.flush();};
            row(progress);if(time.steps%1000==0)std::cout<<"PROGRESS step="<<time.steps<<" stage="<<stage<<" shear="<<next.shear<<" Fs="<<shear_force<<" residual="<<residual<<" KE="<<kr<<" damage="<<stat.max_damage<<" dt="<<dt<<std::endl;
            bool good=stage_time>=ramp&&residual<=o.force_tolerance&&imbalance<=o.balance_tolerance&&kr<=o.kinetic_tolerance;checks=good?checks+1:0;
            if(checks>=3){row(eq);snapshot("equilibrium_"+std::to_string(stage));std::cout<<"ACCEPTED stage="<<stage<<" shear="<<next.shear<<" Fs="<<shear_force<<" damage="<<stat.max_damage<<std::endl;++stage;stage_time=0;checks=0;}
        }
        if(stage_time>box.maximum_stage_duration_s)throw std::runtime_error("Direct-shear stage did not equilibrate; no completed curve claimed");
        if(time.steps%200==0&&(stage!=prior_stage||recovery::checkpoint.due(time.steps,stage!=prior_stage)||std::chrono::steady_clock::now()-lastsave>=std::chrono::seconds(box.checkpoint_wall_seconds))){save();lastsave=std::chrono::steady_clock::now();}
        if(recovery::checkpoint.stop_step&&time.steps>=recovery::checkpoint.stop_step){save();return 86;}
    }
    std::ofstream summary(out/"summary.json");summary<<"{\"status\":\"COMPLETE_LOCAL_DIRECT_SHEAR_FUNCTIONAL_CASE\",\"normal_stages\":"<<box.normal_increments<<",\"shear_stages\":"<<box.shear_increments<<",\"physical_mass_preserved\":true,\"material_parameters_changed\":false,\"bulk_experiment_validation_claimed\":false,\"self_healing_strength_recovery_claimed\":false}\n";
    return 0;
}

int main(int argc,char**argv)
{
    try{std::string cells,sand,boxfile;std::vector<char*>args{argv[0]};for(int i=1;i<argc;++i){std::string a=argv[i];
      if(a.rfind("--snsni-cells=",0)==0)cells=a.substr(14);else if(a.rfind("--sand-input=",0)==0)sand=a.substr(13);else if(a.rfind("--box-config=",0)==0)boxfile=a.substr(13);else if(recovery::checkpoint.take(a)){}else args.push_back(argv[i]);}
      if(cells.empty()||sand.empty()||boxfile.empty())throw std::invalid_argument("cells, sand-input and box-config are required");
      HybridOptions options;options.parseHybrid(int(args.size()),args.data());ShearBoxOptions box;box.read(boxfile);
      const PeriExecution execution(options.execution);
      tbb::task_arena arena(int(options.threads));return arena.execute([&]{return runDirectShear(options,cells,sand,box,execution);});
    }catch(const std::exception&e){std::cerr<<"DIRECT_SHEAR_ERROR "<<e.what()<<std::endl;return 1;}
}
