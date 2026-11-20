/** S-NSNI/RKPM DPC core + layered NOSB-PD shell with optional transport exports.
 * One common explicit displacement/velocity state, physical lumped masses,
 * common-trace vector penalty and frictionless reconstructed outer contact.
 * No WCSPH fluid, material retuning, particle erosion or force rescaling.
 */
#include "pd_snsni_capsule_fixture.h"
#include "solid_snsni_plastic.h"
#include "pd_radial_layered_surface.h"
#include "pd_column_face_surface.h"
#include "recovery_checkpoint.h"
#include "pd_extrinsic_crack_growth.h"
#include "pd_fracture_topology_checkpoint.h"
using namespace SPH::solid_dynamics;
using namespace SPH::peridynamics;
using Field=StdVec<Vecd>;
std::string fracture_bank_input;
bool use_column_interface=false,column_handoff_only=false,local_platen_step_bound=false,fracture_contact_only=false;
bool affine_crack_contact=false,affine_contact_handoff_only=false;
bool affine_inner_interface=false,inner_interface_handoff_only=false;

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
    scalar("BrokenBondFraction",damage);scalar("Pressure_Pa",pressure);scalar("J",J);scalar("PlasticShear",plastic);scalar("FractureDissipationDensity",dissipation);vector("Velocity",v);vector("CrackNormalReference",normal);
    f<<"</PointData></Piece></PolyData></VTKFile>\n";if(!f)throw std::runtime_error("VTP output failed");
}

void vtkPlates(const fs::path&path,Real R,Real compression)
{
    std::ofstream f(path);f<<std::setprecision(17)<<"<VTKFile type=\"PolyData\" version=\"0.1\"><PolyData><Piece NumberOfPoints=\"8\" NumberOfPolys=\"2\"><Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
    for(int side:{-1,1})for(auto xy:{std::pair{-1,-1},std::pair{1,-1},std::pair{1,1},std::pair{-1,1}})f<<xy.first*1.35*R<<' '<<xy.second*1.35*R<<' '<<side*(R-compression/2)<<' ';
    f<<"</DataArray></Points><Polys><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">0 1 2 3 4 5 6 7</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">4 8</DataArray></Polys></Piece></PolyData></VTKFile>\n";
}

int runSnsniCapsule(const HybridOptions&o,const std::string&cells,const PeriExecution&execution)
{
    if(!o.layered_kinematics||o.shell_material!="finite_volume_biot"||o.fracture_model!="extrinsic_cohesive"||o.fracture_kinematics!="physical_cohesive"||o.plate_friction!=0||fracture_bank_input.empty())
        throw std::invalid_argument("This entry requires layered NOSB, unchanged Biot elasticity, extrinsic_cohesive/physical_cohesive fracture and --fracture-facets");
    fs::path out=fs::absolute(o.output);if(fs::exists(out)&&!fs::is_empty(out))throw std::runtime_error("Output must be new or empty; no overwrite");
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
                         execution.makeSpatialDerivatives());
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

    PDLayeredFractureReference sr(relation,geometry);
    InteractionWithUpdate<NosbPDQuadraticSecondStep> quadratic_second(relation,sr);
    InteractionDynamics<NosbPDQuadraticThirdStep> quadratic_third(relation,sr);
    PDBrokenBondEnergy sb(relation,sr,o.shell_E,o.shell_nu,o.pd_hourglass,shell.getSPHAdaptation().getKernel()->CutOffRadius(),
                          execution.makeSpatialDerivatives());
    system.initializeSystemCellLinkedLists();system.initializeSystemConfigurations();shape.exec();sr.initialize();sb.rebuild();
    SolidForcePort core_port(core,SolidForcePort::Storage::ForcePrior),pd_port(shell,SolidForcePort::Storage::PDAccelerationPrior);
    StdVec<Vecd>points,outer_points;StdVec<Real>areas;
    for(const auto&face:cr.faces)if(face.boundary)for(const auto&q:face.q){points.push_back(q.x);areas.push_back(q.w);}
    for(const auto&n:sphereDirections(R,o.contact_quadrature_scale*o.spacing))outer_points.push_back(R*n);
    // The reference boundary coordinate is the material attachment. A nearest
    // node is only its initial search anchor and may lie across a later crack.
    // Resolve the same point's visible side after a cut; never average fragments.
    using Cell=std::array<long long,3>;struct Hash{size_t operator()(const Cell&a)const{size_t h=0;for(auto v:a)h^=std::hash<long long>{}(v)+0x9e3779b9+(h<<6)+(h>>2);return h;}};
    std::unordered_map<Cell,StdVec<size_t>,Hash>grid;Real cell_size=2*o.pd_spacing;
    auto cell=[&](const Vecd&x){return Cell{(long long)std::floor(x[0]/cell_size),(long long)std::floor(x[1]/cell_size),(long long)std::floor(x[2]/cell_size)};};
    for(size_t j=0;j<sn;++j)grid[cell(sX[j])].push_back(j);
    auto anchors=[&](const StdVec<Vecd>&targets){StdVec<size_t>ids;for(const auto&q:targets){auto c=cell(q);Real best=std::numeric_limits<Real>::infinity();size_t owner=sn;
        for(int a=-1;a<=1;++a)for(int b=-1;b<=1;++b)for(int d=-1;d<=1;++d){auto it=grid.find(Cell{c[0]+a,c[1]+b,c[2]+d});if(it!=grid.end())for(size_t j:it->second){Real r=(sX[j]-q).squaredNorm();if(r<best){best=r;owner=j;}}}
        if(owner==sn)throw std::runtime_error("No shell material owner for interface/plate trace");ids.push_back(owner);}return ids;};
    Field inner_affine_positions;
    const auto inner_anchors=anchors(points);
    auto makeInnerSurface=[&](bool affine,const Field&target,const Field&velocity){
        using Ownership=VisibleCrackSurface::AnchorOwnership;
        using Transfer=VisibleCrackSurface::TransferOrder;
        return std::make_unique<VisibleCrackSurface>(pd_port,points,inner_anchors,sr,o.pd_spacing,o.interface_shell_support,
            affine?1:2,Field{},target,velocity,Ownership::ReferencePoint,
            affine?Transfer::AllowAffinePositionTransfer:Transfer::PreservePolynomial,affine?&inner_affine_positions:nullptr);
    };
    auto si=makeInnerSurface(affine_inner_interface&&recovery::checkpoint.resume.empty(),{},{});
    VisibleCrackSurface so(pd_port,outer_points,anchors(outer_points),sr,o.pd_spacing,o.interface_shell_support,2,
                           {},{},{},VisibleCrackSurface::AnchorOwnership::ReferencePoint,
                           VisibleCrackSurface::TransferOrder::AllowAffinePositionTransfer);
    Real Eref=o.constraint_E>0?o.constraint_E:o.core_E,Ktie=o.tie_factor*Eref/o.spacing,Kcontact=o.contact_factor*Eref/o.spacing;
    auto cs=std::make_unique<SnsniCoupledSurface>(core_port,points,cr);
    auto tie=std::make_unique<BondedSolidInterface>(*si,*cs,areas,Ktie,R,execution.makeBondedInterface());
    auto installAffineInner=[&](const Field&targets,const Field&velocities){
        auto next=makeInnerSurface(true,targets,velocities);
        tie=std::make_unique<BondedSolidInterface>(*next,*cs,areas,Ktie,R,execution.makeBondedInterface());
        si=std::move(next);
    };
    StdVec<Real>outer_areas(outer_points.size(),4*Pi*R*R/outer_points.size());
    std::unique_ptr<ColumnFaceSolidSurface> column_surface;
    std::unique_ptr<SnsniCoupledSurface> column_core;
    StdVec<Real> column_areas;
    if(use_column_interface){
        column_surface=std::make_unique<ColumnFaceSolidSurface>(pd_port,geometry,-1);
        column_core=std::make_unique<SnsniCoupledSurface>(core_port,column_surface->points,cr);
        const Real total_area=std::accumulate(areas.begin(),areas.end(),Real(0));
        const Real shell_volume=std::accumulate(svol,svol+sn,Real(0));
        // Radial columns of identical thickness have area proportional to volume.
        // Integrate the same total interface area as the original core boundary.
        for(const auto&c:geometry.columns){Real v=0;for(size_t j:c)v+=svol[j];column_areas.push_back(total_area*v/shell_volume);}
    }
    auto installColumnInterface=[&](){
        tie=std::make_unique<BondedSolidInterface>(*column_surface,*column_core,column_areas,Ktie,R,execution.makeBondedInterface());
        si.reset();cs.reset();
    };
    if(use_column_interface&&recovery::checkpoint.resume.empty())installColumnInterface();
    SolidSurfacePlaneContact contact(so,outer_areas,Kcontact,R);
    auto bank=ExtrinsicCrackGrowth::readBank(fracture_bank_input,geometry);
    using SurfaceLaw=ExtrinsicCrackGrowth::SurfaceLaw;
    auto crack_projection=execution.makeCrackProjection();
    ExtrinsicCrackGrowth growth(sr,sb,pd_port,o.fracture_stress,o.fracture_energy,o.pd_spacing,Kcontact,bank,
                                fracture_contact_only?SurfaceLaw::ContactOnly:SurfaceLaw::Cohesive,crack_projection.get());
    growth.affine_contact_traces=affine_crack_contact&&recovery::checkpoint.resume.empty();
    cr.releasePreprocessingCache();
    size_t count=cn+sn;StdVec<Real>mass(count);Field X(count);
    for(size_t j=0;j<cn;++j){mass[j]=cm[j];X[j]=cX[j];}for(size_t j=0;j<sn;++j){mass[cn+j]=sm[j];X[cn+j]=sX[j];}
    SnsniCentralDifference time(mass);Field bulk_force(count,Vecd::Zero()),internal(count,Vecd::Zero()),cohesive_force(sn,Vecd::Zero());
    Field inverse_shell_mass(sn);for(size_t j=0;j<sn;++j)inverse_shell_mass[j]=Vecd::Constant(1/sm[j]);
    auto*sF=sp.getVariableDataByName<Matd>("DeformationGradient");auto*sstress=sp.getVariableDataByName<Matd>("PDStress");
    Real shell_mass=std::accumulate(sm,sm+sn,Real(0)),exact_mass=o.shell_rho*4*Pi*(R*R*R-rc*rc*rc)/3;
    if(std::abs(shell_mass/exact_mass-1)>1e-10)throw std::runtime_error("Physical shell mass mismatch");
    Real radial_spacing=o.thickness;for(size_t k=1;k<rule.xi.size();++k)radial_spacing=std::min(radial_spacing,.5*o.thickness*(rule.xi[k]-rule.xi[k-1]));
    auto&mat=DynamicCast<HughesWingetSolid>(&shell,shell.getMatterMaterial());Real radial_dt=o.pd_cfl*std::min(radial_spacing,o.pd_spacing)/mat.ReferenceSoundSpeed(),core_dt=o.cfl/std::sqrt(cb.linear_frequency_bound_squared);
    auto setPositions=[&](const Field&d,const Field&v){for(size_t j=0;j<cn;++j){cx[j]=cX[j]+d[j];cv[j]=v[j];}for(size_t j=0;j<sn;++j){sx[j]=sX[j]+d[cn+j];sv[j]=v[cn+j];}};
    auto assembleShellAndConstraints=[&](Real compression){
        sb.evaluateAll();sb.recomputeForces();
        core_port.clear();pd_port.clear();tie->apply();contact.apply(-R+compression/2,R-compression/2);growth.applyContact();
        for(size_t j=0;j<cn;++j){internal[j]=cb.internal[j];bulk_force[j]=core_port.external[j]-internal[j];}
        for(size_t j=0;j<sn;++j){internal[cn+j]=-sb.internal[j]-sb.hourglass[j];bulk_force[cn+j]=sm[j]*pd_port.external[j]-internal[cn+j];}
    };
    auto forceAt=[&](const Field&d,const Field&v,Real compression,bool commit){setPositions(d,v);Field cd(d.begin(),d.begin()+cn);cb.evaluate(cd,o.core_rho,o.cfl);if(commit)cb.commit();else cb.discard();assembleShellAndConstraints(compression);};
    auto recoverable=[&](){return sb.energy()+growth.stored+growth.contact_energy+tie->energy+contact.energy;};
    size_t stage=1,checks=0;Real stage_time=0,elapsed=0,Wp=0,Wcore=0,Wd=0,Walgorithmic=0,interface_switch_work=0,inner_hold_elapsed=0;bool onset=false,fully_failed=false,column_hold=false,inner_hold=false;
    auto compression=[&](Real t){Real a=column_hold?0:std::clamp((inner_hold?stage_time:t)/o.ramp_duration,Real(0),Real(1));return o.diameter*o.strain*(Real(stage-1)+.5*(1-std::cos(Pi*a)))/o.load_steps;};
    forceAt(time.displacement,time.velocity,0,false);Real U0=recoverable(),K0=time.kinetic();
    std::ofstream parameters(out/"parameters.json");parameters<<std::setprecision(17)<<"{\"core_backend\":\"SNSNI_RKPM_DPC_unchanged\",\"shell_backend\":\"layered_NOSB_PD_live_reconstruction\",\"shell_constitutive\":\"unchanged_finite_volume_Biot\",\"fracture_representation\":"<<std::quoted(fracture_contact_only?"strength_triggered_bond_deletion_contact_only":"strength_triggered_extrinsic_cohesive_true_bonds")<<",\"cohesive_integrator\":"<<std::quoted(fracture_contact_only?"disabled":"explicit_bulk_half_kicks_implicit_objective_surface_drift")<<",\"Gc_used_for_surface_softening\":"<<(fracture_contact_only?"false":"true")<<",\"mixed_mode_assumption\":"<<std::quoted(fracture_contact_only?"equal_strength":"equal_strength_and_Gc")<<",\"candidate_geometry\":\"through_thickness_radial_Voronoi_facets\",\"E_Pa\":"<<o.shell_E<<",\"nu\":"<<o.shell_nu<<",\"ft_Pa\":"<<o.fracture_stress<<",\"Gc_J_m2\":"<<o.fracture_energy<<",\"delta_f_m\":"<<growth.delta_f<<",\"core_particles\":"<<cn<<",\"shell_particles\":"<<sn<<",\"candidate_facets\":"<<bank.size()<<",\"shell_mass_kg\":"<<shell_mass<<",\"interface_K_Pa_per_m\":"<<Ktie<<",\"platen_and_crack_contact_K_Pa_per_m\":"<<Kcontact<<",\"execution_backend\":"<<std::quoted(execution.backendName())<<",\"crack_projection_backend\":"<<std::quoted(crack_projection?"cuda_fp64":"cpu_fp64")<<",\"crack_local_contact_backend\":"<<std::quoted(crack_projection&&crack_projection->assemblesLocalContact()?"cuda_fp64":"cpu_fp64")<<",\"crack_node_assembly_backend\":"<<std::quoted(crack_projection&&crack_projection->assemblesTracePullback()?"cuda_fp64":"cpu_fp64")<<",\"shell_derivatives_backend\":"<<std::quoted(sb.hasDeviceDerivatives()?"cuda_fp64":"cpu_fp64")<<",\"core_derivatives_backend\":"<<std::quoted(cb.hasDeviceDerivatives()?"cuda_fp64":"cpu_fp64")<<",\"shell_stress_force_backend\":"<<std::quoted(sb.hasDeviceDerivatives()?"cuda_fp64":"cpu_fp64")<<",\"shell_stabilization_backend\":"<<std::quoted(sb.hasDeviceDerivatives()?"cuda_fp64":"cpu_fp64")<<",\"core_internal_force_backend\":"<<std::quoted(cb.hasDeviceDerivatives()?"cuda_fp64":"cpu_fp64")<<",\"bonded_interface_backend\":"<<std::quoted(execution.bondedInterfaceOnDevice()?"cuda_fp64":"cpu_fp64")<<",\"core_parameters_changed\":false,\"source_hash\":"<<std::quoted(CAPSULE_SOURCE_SHA256)<<"}\n";parameters.close();
    if(o.preflight_only){Field d(count),v(count);Matd A=Matd::Identity();A(0,0)=.998;A(1,1)=1.001;A(2,2)=.999;for(size_t j=0;j<count;++j){d[j]=(A-Matd::Identity())*X[j];if(j<cn)d[j]+=Vecd(5e-8*rc/.0005,0,0);v[j]=Vecd(.2*std::sin(37*X[j][0]/R),-.1,.3*X[j][2]/R);}forceAt(d,v,.004*2*R,false);
        Real affine=0,native=0;for(size_t l=0;l<cr.cells.size();++l){affine=std::max(affine,(cb.F[l]-A).norm());auto a=cb.model.update(cb.F[l],cb.history[l]);native=std::max(native,(a.sigma-cb.pending[l].material.sigma).norm());}for(size_t j=0;j<sn;++j)affine=std::max(affine,(sF[j]-A).norm());
        Real power=0;for(size_t j=0;j<cn;++j)power+=core_port.external[j].dot(v[j]);for(size_t j=0;j<sn;++j)power+=sm[j]*pd_port.external[j].dot(v[cn+j]);Real eps=1e-12;Field plus=d,minus=d;for(size_t j=0;j<count;++j){plus[j]+=eps*v[j];minus[j]-=eps*v[j];}forceAt(plus,v,.004*2*R,false);Real ep=tie->energy+contact.energy;forceAt(minus,v,.004*2*R,false);Real em=tie->energy+contact.energy;Real derivative=(ep-em)/(2*eps),error=std::abs(derivative+power)/std::max({Real(1e-20),std::abs(derivative),std::abs(power)});bool pass=affine<1e-8&&native==0&&error<1e-5;
        std::ofstream f(out/"preflight.json");f<<std::setprecision(17)<<"{\"status\":\""<<(pass?"PASSED":"FAILED")<<"\",\"affine_F_error\":"<<affine<<",\"native_DPC_stress_difference_Pa\":"<<native<<",\"constraint_energy_gradient_relative_error\":"<<error<<",\"fracture_contact_only\":"<<(fracture_contact_only?"true":"false")<<",\"affine_crack_contact\":"<<(affine_crack_contact?"true":"false")<<",\"affine_inner_interface\":"<<(affine_inner_interface?"true":"false")<<",\"physical_histories_advanced\":0}\n";return pass?0:2;}
    std::ofstream progress(out/"progress.csv"),eq(out/"equilibria.csv"),events(out/"fracture_events.csv");
    const char*header="step,stage,elapsed_s,strain,compression_m,dt_s,upper_N,lower_N,force_residual,plate_imbalance,kinetic_ratio,core_min_J,shell_min_J,broken_bond_pairs,crack_facets,failed_cohesive_quadrature,max_crack_gap_m,cohesive_stored_J,fracture_work_J,shell_stored_J,interface_energy_J,contact_energy_J,crack_contact_energy_J,max_penetration_m,platen_work_J,core_internal_work_J,damping_loss_J,cohesive_algorithmic_loss_J,stabilizer_transfer_loss_J,insertion_defect_J,mechanical_budget_J\n";
    progress<<std::setprecision(17)<<header;eq<<std::setprecision(17)<<header;events<<std::setprecision(17)<<"event,step,stage,strain,candidate,strength_ratio,broken_bond_pairs,max_gap_m\n";
    std::ofstream surface_audit(out/"surface_reconstruction_events.csv");surface_audit<<std::setprecision(17)<<"step,stage,strain,candidate,inner_reassigned_points,outer_reassigned_points,maximum_affine_error,maximum_polynomial_error,maximum_trace_l1,insertion_energy_fraction,core_dt_s,shell_dt_s,stabilizer_dt_s,tie_dt_s,platen_dt_s,crack_contact_dt_s,platen_position_fallback_points,platen_expanded_support_points\n";
std::vector<std::pair<Real,std::string>>frame_index;
    auto snapshot=[&](const std::string&tag,Real u){Field cpos(cr.cells.size()),cd(time.displacement.begin(),time.displacement.begin()+cn);StdVec<Real>pressure(cr.cells.size()),J(cr.cells.size()),alpha(cr.cells.size());
        for(size_t l=0;l<cr.cells.size();++l){cpos[l]=cr.cells[l].center+cr.value(l,cd);pressure[l]=cb.pending[l].material.p;J[l]=cb.F[l].determinant();alpha[l]=cb.history[l].alpha;}vtkPoints(out/"vtp"/("core_"+tag+".vtp"),cpos,{}, {},pressure,J,alpha,{},{});
        StdVec<Real>broken(sn),shell_pressure(sn),shell_J(sn);for(size_t j=0;j<sn;++j){size_t dead=0;for(size_t n=0;n<sr.rows[j].size();++n)dead+=!sr.live(j,n);broken[j]=Real(dead)/sr.rows[j].size();shell_pressure[j]=-sstress[j].trace()/3;shell_J[j]=sF[j].determinant();}
        vtkPoints(out/"vtp"/("shell_"+tag+".vtp"),Field(sx,sx+sn),Field(sv,sv+sn),broken,shell_pressure,shell_J,{}, {},{});vtkPlates(out/"vtp"/("platens_"+tag+".vtp"),R,u);
        std::ofstream faces(out/"vtp"/("cracks_"+tag+".vtp"));faces<<std::setprecision(17)<<"<VTKFile type=\"PolyData\" version=\"0.1\"><PolyData><Piece NumberOfPoints=\""<<8*growth.patches.size()<<"\" NumberOfPolys=\""<<2*growth.patches.size()<<"\"><Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
        for(const auto&p:growth.patches)for(auto*s:{p.minus.get(),p.plus.get()})for(size_t k=0;k<4;++k)faces<<s->position(k).transpose()<<' ';faces<<"</DataArray></Points><Polys><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">";for(size_t k=0;k<2*growth.patches.size();++k)faces<<4*k<<' '<<4*k+1<<' '<<4*k+2<<' '<<4*k+3<<' ';faces<<"</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">";for(size_t k=0;k<2*growth.patches.size();++k)faces<<4*(k+1)<<' ';faces<<"</DataArray></Polys><CellData>";
        // Two geometrical faces per facet. Bond deletion, spent Gc and actual
        // opening are different observables; exporting them separately avoids
        // treating a point damage colour as evidence of a separated crack.
        StdVec<Real>work_fraction,failed_fraction,normal_opening;
        if(fracture_contact_only){
            for(const auto&patch:growth.patches){
                failed_fraction.push_back(1);
                normal_opening.push_back(growth.maximumNormalOpening(patch));
            }
        }else{
        size_t history_index=0;auto shell_positions=growth.positions();CohesiveSurfaceEvolution fracture_view(growth.ft,growth.Gc,growth.quadrature());
        for(const auto&p:growth.patches){const auto&candidate=bank[p.candidate];Real area=0,spent=0,failed=0,gap=0;for(size_t q=0;q<candidate.point.size();++q,++history_index){Real a=candidate.area[q],k=growth.history[history_index].kappa;area+=a;spent+=a*std::min(k/growth.delta_f,Real(1));if(k>=growth.delta_f)failed+=a;gap=std::max(gap,std::max(fracture_view.geometry(history_index,shell_positions).normal_gap,Real(0)));}work_fraction.push_back(spent/area);failed_fraction.push_back(failed/area);normal_opening.push_back(gap);}
        }
        auto cell_data=[&](const char*name,const StdVec<Real>&values){faces<<"<DataArray type=\"Float64\" Name=\""<<name<<"\" format=\"ascii\">";for(Real value:values)faces<<value<<' '<<value<<' ';faces<<"</DataArray>";};
        if(!fracture_contact_only)cell_data("FractureWorkFraction",work_fraction);cell_data("FullyFailedAreaFraction",failed_fraction);cell_data("MaximumNormalOpening_m",normal_opening);faces<<"</CellData></Piece></PolyData></VTKFile>\n";
        std::ofstream vtm(out/"vtp"/("capsule_"+tag+".vtm"));vtm<<"<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\"><vtkMultiBlockDataSet>";size_t part=0;for(auto name:{"shell","core","platens","cracks"})vtm<<"<DataSet index=\""<<part++<<"\" name=\""<<name<<"\" file=\""<<name<<'_'<<tag<<".vtp\"/>";vtm<<"</vtkMultiBlockDataSet></VTKFile>\n";
        faces.close();vtm.close();if(!faces||!vtm)throw std::runtime_error("Crack trajectory publication failed");
        frame_index.push_back({elapsed,tag});std::ofstream index(out/"capsule.pvd.tmp");index<<std::setprecision(17)<<"<VTKFile type=\"Collection\" version=\"0.1\"><Collection>\n";for(const auto&f:frame_index)index<<"<DataSet timestep=\""<<f.first<<"\" file=\"vtp/capsule_"<<f.second<<".vtm\"/>\n";index<<"</Collection></VTKFile>\n";index.close();if(!index)throw std::runtime_error("Trajectory index publication failed");fs::rename(out/"capsule.pvd.tmp",out/"capsule.pvd");
    };
    if(recovery::checkpoint.resume.empty())snapshot("0",0);
    auto refresh=[&](){if(si){si->rebuildChanged();tie->refreshTimeStepBound();}so.rebuildChanged();contact.refreshTimeStepBound();};
    StdVec<Vecd>old_si_transfer,old_so_transfer,old_si_velocity,old_so_velocity;
    Field previous_inner_affine_positions;
    bool surface_transfer_prepared=false;
    auto prepareTransfer=[&](){
        if(si){old_si_transfer=si->transferValues();old_si_velocity=si->transferVelocities();si->captureTransferValues();}
        old_so_transfer=so.transferValues();old_so_velocity=so.transferVelocities();so.captureTransferValues();
        previous_inner_affine_positions=std::move(inner_affine_positions);
        surface_transfer_prepared=true;
    };
    auto rollbackTransfer=[&](){
        if(!surface_transfer_prepared)return;
        inner_affine_positions=std::move(previous_inner_affine_positions);
        if(si)si->restoreTransferValues(old_si_transfer,old_si_velocity);
        so.restoreTransferValues(old_so_transfer,old_so_velocity);
    };
    auto eventEnergy=[&](){assembleShellAndConstraints(compression(stage_time));growth.statistics(growth.positions());return recoverable();};
    using CrackArchive=FractureTopologyCheckpoint;std::vector<uint64_t>point_ids(sn);for(size_t j=0;j<sn;++j)point_ids[j]=j;
    auto transfer=[&](recovery::Archive&a){
        const auto version=a.reading?recovery::checkpoint.read_version:recovery::checkpoint.schema_version;
        if(version==5||version==7||version==9){
            std::string model="fixed_column_trace";a.fields(model,interface_switch_work,column_hold);
            if(model!="fixed_column_trace"||!use_column_interface||!std::isfinite(interface_switch_work))
                throw std::runtime_error("Invalid fixed-column checkpoint interface");
            if(a.reading)installColumnInterface();
        }
        a.fields(stage,stage_time,elapsed,checks,onset,fully_failed,Wp,Wcore,Wd,Walgorithmic,U0,K0,bulk_force,internal,cohesive_force);
        recovery::integrator(a,time);recovery::solid(a,cb);recovery::contact(a,contact);a.fields(growth.insertion_defect,growth.max_insertion_fraction,growth.initial_dissipation,growth.insertion_events,growth.rejected_events,growth.stabilizer_transfer_loss,growth.unresolved_candidates);
        StdVec<Vecd> inner_transfer,inner_velocity;
        if(si){inner_transfer=si->transferValues();inner_velocity=si->transferVelocities();}
        StdVec<Vecd> outer_transfer=so.transferValues(),outer_velocity=so.transferVelocities();
        StdVec<StdVec<Vecd>> minus_transfer,plus_transfer,minus_velocity,plus_velocity;
        if(!a.reading)for(const auto&p:growth.patches){minus_transfer.push_back(p.minus->transferValues());plus_transfer.push_back(p.plus->transferValues());minus_velocity.push_back(p.minus->transferVelocities());plus_velocity.push_back(p.plus->transferVelocities());}
        a.fields(inner_transfer,outer_transfer,minus_transfer,plus_transfer,inner_velocity,outer_velocity,minus_velocity,plus_velocity);
        std::vector<size_t>order;CrackArchive::Bytes bytes;auto traction=growth.last_traction;if(!a.reading)for(const auto&p:growth.patches)order.push_back(p.candidate);a.value(order);a.value(traction);auto candidate_flags=growth.active;a.value(candidate_flags);
        CrackArchive::Histories history;CrackArchive::IMEXBlock block;block.name="shell_extrinsic_equal_mode";block.ft=growth.ft;block.Gc=growth.Gc;size_t qindex=0;
        for(size_t k:order){if(k>=bank.size())throw std::runtime_error("Checkpoint candidate ID outside source bank");for(size_t q=0;q<bank[k].point.size();++q){block.history_ids.push_back(64*bank[k].id+q);block.history.push_back(a.reading?CohesiveInterfaceIMEX::History{}:growth.history.at(qindex));++qindex;}}history.imex.push_back(block);
        if(!a.reading)bytes=CrackArchive::encode(sr,recovery::checkpoint.fingerprint,point_ids,history);a.value(bytes);
        if(a.reading&&fracture_contact_only&&version<6&&(!order.empty()||onset||fully_failed))
            throw std::runtime_error("The change to contact-only fracture requires an intact checkpoint before the first cut");
        // The frozen numerical operator is history too. Store its exact rows,
        // including the intact initial operator, instead of relying on algebraic
        // reassociation during reconstruction to reproduce zero forces.
        StdVec<StdVec<Vecd>>saved_gradient;if(!a.reading){saved_gradient.resize(sn);for(size_t i=0;i<sn;++i)for(const auto&e:sr.rows[i])saved_gradient[i].push_back(e.gradient);}Real saved_elastic_dt=sr.elastic_dt_bound;a.fields(saved_gradient,saved_elastic_dt);
        a.span(cx,cn);a.span(cv,cn);a.span(sx,sn);a.span(sv,sn);
        if(version>=8){
            a.value(growth.affine_contact_transfer_positions);
            const auto&positions=growth.affine_contact_transfer_positions;
            if(!affine_crack_contact||(!positions.empty()&&positions.size()!=sn)||
                std::any_of(positions.begin(),positions.end(),[](const Vecd&q){return !q.allFinite();}))
                throw std::runtime_error("Invalid affine crack contact checkpoint transfer state");
            growth.affine_contact_traces=true;
        }
        if(version>=10){
            a.fields(inner_affine_positions,inner_hold,inner_hold_elapsed);
            if(!std::isfinite(inner_hold_elapsed)||inner_hold_elapsed<0)
                throw std::runtime_error("Invalid affine inner-interface hold clock");
            if(!affine_inner_interface||(!inner_affine_positions.empty()&&inner_affine_positions.size()!=sn)||
                std::any_of(inner_affine_positions.begin(),inner_affine_positions.end(),[](const Vecd&q){return !q.allFinite();}))
                throw std::runtime_error("Invalid affine core-shell checkpoint transfer state");
        }
        if(a.reading){CrackArchive::load(sr,bytes,recovery::checkpoint.fingerprint,point_ids,history);
            if(fracture_contact_only&&version<6&&sr.physical_broken_bonds!=0)
                throw std::runtime_error("The original checkpoint already contains broken bonds");
            if(saved_gradient.size()!=sn||!std::isfinite(saved_elastic_dt)||saved_elastic_dt<=0)throw std::runtime_error("Checkpoint live-operator dimensions/bound");
            for(size_t i=0;i<sn;++i){if(saved_gradient[i].size()!=sr.rows[i].size())throw std::runtime_error("Checkpoint live-operator row count");Matd affine=Matd::Zero();for(size_t k=0;k<sr.rows[i].size();++k){const auto&g=saved_gradient[i][k];if(!g.allFinite()||(!sr.live(i,k)&&g.squaredNorm()!=0))throw std::runtime_error("Invalid checkpoint live gradient");sr.rows[i][k].gradient=g;affine+=(sX[sr.rows[i][k].j]-sX[i])*g.transpose();}if((affine-Matd::Identity()).norm()>1e-8)throw std::runtime_error("Checkpoint live operator lost affine completeness");}sr.elastic_dt_bound=saved_elastic_dt;
            if(version>=10)installAffineInner(inner_transfer,inner_velocity);
            else if(si)si->restoreTransferValues(inner_transfer,inner_velocity);
            else if(!inner_transfer.empty()||!inner_velocity.empty())throw std::runtime_error("Fixed-column restart includes a dynamic inner trace");
            so.restoreTransferValues(outer_transfer,outer_velocity);growth.restorePatches(order,history.imex[0].history,minus_transfer,plus_transfer,minus_velocity,plus_velocity);if(traction.size()!=growth.history.size())throw std::runtime_error("Checkpoint adhesion-traction dimensions");growth.last_traction=std::move(traction);if(candidate_flags.size()!=growth.active.size())throw std::runtime_error("Checkpoint candidate-state dimensions");for(size_t k=0;k<candidate_flags.size();++k){if(candidate_flags[k]>2||(candidate_flags[k]==1)!=(growth.active[k]==1))throw std::runtime_error("Checkpoint candidate-state mismatch");}growth.active=std::move(candidate_flags);sb.rebuild();refresh();}
        if(time.displacement.size()!=count||cb.history.size()!=cr.cells.size())throw std::runtime_error("Checkpoint geometry mismatch");
    };
    auto save=[&](){progress.flush();eq.flush();events.flush();recovery::checkpoint.save(out,{"progress.csv","equilibria.csv","fracture_events.csv"},time.steps,transfer);};
    StdVec<size_t> bounded_rows;
    Real platen_bound=std::numeric_limits<Real>::infinity();
    size_t bound_version=~size_t(0);
    auto boundRows=[&](const Field* trial,Real trial_u){
        StdVec<size_t> rows;
        const Real top=R-.5*trial_u,bottom=-top;
        StdVec<Real> trial_z;
        if(trial){
            trial_z.resize(sn);
            for(size_t j=0;j<sn;++j)trial_z[j]=sX[j][2]+(*trial)[cn+j][2];
        }
        for(size_t i=0;i<so.rows.size();++i){
            const auto&r=so.rows[i];
            if(r.index.empty())continue;
            const Real old_z=contact.position[i][2];
            Real z=old_z,margin=0;
            if(trial){
                z=0;
                Real magnitude=std::abs(old_z)+std::abs(top);
                for(size_t a=0;a<r.index.size();++a){
                    const Real value=r.weight[a]*trial_z[r.index[a]];
                    z+=value;
                    magnitude+=std::abs(value);
                }
                margin=(2*r.index.size()+4)*std::numeric_limits<Real>::epsilon()*magnitude;
            }
            if(!std::isfinite(z)||!std::isfinite(margin))throw std::runtime_error("Invalid platen trial geometry");
            if(std::min(old_z,z)-margin<=bottom||std::max(old_z,z)+margin>=top)rows.push_back(i);
        }
        if(bound_version!=so.mappingVersion()||rows!=bounded_rows){
            platen_bound=ReferenceSolidSurface::stableBound(so,nullptr,contact.stiffness,.2,&rows);
            bounded_rows=std::move(rows);bound_version=so.mappingVersion();
        }
        return platen_bound;
    };
    if(!recovery::checkpoint.resume.empty()){progress.close();eq.close();events.close();recovery::checkpoint.load(out,{"progress.csv","equilibria.csv","fracture_events.csv"},transfer);progress.open(out/"progress.csv",std::ios::app);eq.open(out/"equilibria.csv",std::ios::app);events.open(out/"fracture_events.csv",std::ios::app);progress<<std::setprecision(17);eq<<std::setprecision(17);events<<std::setprecision(17);
        Field expected=bulk_force;assembleShellAndConstraints(compression(stage_time));Real error=0,scale=0;for(size_t j=0;j<count;++j){error=std::max(error,(bulk_force[j]-expected[j]).norm());scale=std::max(scale,expected[j].norm());}if(error>1e-9*std::max(scale,Real(1e-20)))throw std::runtime_error("Resumed physical-fracture force does not reproduce committed checkpoint");
        if(use_column_interface&&si){
            if(stage_time!=0||stage<2||checks!=0)
                throw std::runtime_error("The fixed-column interface must be initialized at an accepted load point");
            const Real u=compression(stage_time),old_energy=tie->energy,old_dt=tie->stable_dt;
            const Real old_reaction=.5*(contact.upper+contact.lower);
            installColumnInterface();assembleShellAndConstraints(u);
            interface_switch_work=tie->energy-old_energy;
            const Real energy_fraction=std::abs(interface_switch_work)/std::max(std::abs(Wp),Real(1e-30));
            const Real reaction=.5*(contact.upper+contact.lower);
            Real imbalance=0,redistribution=0;
            for(size_t j=0;j<count;++j){
                Vecd f=bulk_force[j];if(j>=cn)f+=cohesive_force[j-cn];
                imbalance+=f.squaredNorm();redistribution+=(bulk_force[j]-expected[j]).squaredNorm();
            }
            imbalance=std::sqrt(imbalance)/std::max(reaction,Real(1e-30));
            const bool passed=std::isfinite(imbalance)&&energy_fraction<.005&&
                std::abs(reaction-old_reaction)<=1e-10*std::max(old_reaction,Real(1e-30));
            // The old state is still the starting state. A change of nodal
            // load distribution can require relaxation at the same load.
            column_hold=passed&&imbalance>o.force_tolerance;
            std::ofstream audit(out/"column_interface_handoff.json");audit<<std::setprecision(17)
                <<"{\"status\":\""<<(passed?"PASSED":"FAILED")<<"\",\"step\":"<<time.steps<<",\"strain\":"<<u/o.diameter
                <<",\"old_tie_energy_J\":"<<old_energy<<",\"new_tie_energy_J\":"<<tie->energy
                <<",\"discrete_change_work_J\":"<<interface_switch_work<<",\"energy_change_over_external_work\":"<<energy_fraction
                <<",\"old_tie_dt_s\":"<<old_dt<<",\"new_tie_dt_s\":"<<tie->stable_dt
                <<",\"old_reaction_N\":"<<old_reaction<<",\"new_reaction_N\":"<<reaction
                <<",\"nodal_force_change_l2_N\":"<<std::sqrt(redistribution)<<",\"force_residual\":"<<imbalance
                <<",\"force_tolerance\":"<<o.force_tolerance<<",\"virtual_points\":"<<column_areas.size()
                <<",\"hold_at_original_displacement\":"<<(column_hold?"true":"false")<<"}\n";
            if(!passed)throw std::runtime_error("Fixed-column interface handoff exceeds the reaction or work acceptance bound");
            if(column_hold)std::cout<<"FIXED_COLUMN_REEQUILIBRATION_STARTED strain="<<u/o.diameter<<" residual="<<imbalance<<std::endl;
        }
        if(affine_crack_contact&&recovery::checkpoint.read_version<8){
            if(!growth.contact_network)throw std::runtime_error("Affine contact handoff requires existing crack faces");
            const Real u=compression(stage_time),reaction=.5*(contact.upper+contact.lower);
            const Real platen_dt=local_platen_step_bound?boundRows(nullptr,u):contact.stable_dt;
            auto combinedStep=[&](Real crack_dt){
                const Real shell_dt=sb.finiteStrainTimeStep();
                const Real rate=1/(core_dt*core_dt)+1/(shell_dt*shell_dt)+1/(sb.stabilizer_dt*sb.stabilizer_dt)
                    +1/(tie->stable_dt*tie->stable_dt)+1/(platen_dt*platen_dt)+1/(crack_dt*crack_dt);
                return std::min({o.max_dt,1/std::sqrt(rate),cb.acoustic_dt,radial_dt});
            };
            const Real old_force=growth.contact_network->normal_force_sum,old_energy=growth.contact_energy;
            const Real old_gap=growth.max_gap,old_penetration=growth.max_penetration,old_dt=growth.contact_dt;
            const Real old_step=combinedStep(old_dt);
            const size_t old_active=growth.contact_network->active_contacts,old_commits=cb.commits;
            Field old_surface_positions;
            StdVec<Real>old_weights,old_openings;
            auto weightNorm=[](const BilinearCrackSurface&side){
                Real maximum=0;
                for(size_t q=0;q<side.natural.size();++q){
                    Real sum=0;
                    for(Real w:side.rows[4+q].weight)sum+=std::abs(w);
                    maximum=std::max(maximum,sum);
                }
                return maximum;
            };
            for(const auto&patch:growth.patches){
                old_openings.push_back(growth.maximumNormalOpening(patch));
                for(auto*face:{patch.minus.get(),patch.plus.get()}){
                    old_weights.push_back(weightNorm(*face));
                    for(size_t q=0;q<4;++q)old_surface_positions.push_back(face->position(q));
                }
            }
            growth.transferContactTracesToAffine();
            assembleShellAndConstraints(u);
            const auto&network=*growth.contact_network;
            Real corner_error=0,opening_error=0,new_weight_max=0;
            size_t i=0,surface_index=0,patch_index=0;
            std::ofstream weights(out/"affine_crack_trace_weights.csv");
            weights<<std::setprecision(17)<<"candidate,side,old_weight_l1,new_weight_l1\n";
            for(const auto&patch:growth.patches){
                opening_error=std::max(opening_error,std::abs(growth.maximumNormalOpening(patch)-old_openings[patch_index++]));
                for(auto*face:{patch.minus.get(),patch.plus.get()}){
                    const Real new_weight=weightNorm(*face);
                    new_weight_max=std::max(new_weight_max,new_weight);
                    weights<<bank[patch.candidate].id<<','<<(face==patch.minus.get()?"minus":"plus")<<','<<old_weights[surface_index++]<<','<<new_weight<<'\n';
                    for(size_t q=0;q<4;++q)corner_error=std::max(corner_error,(face->position(q)-old_surface_positions[i++]).norm());
                }
            }
            const Real new_step=combinedStep(growth.contact_dt);
            const Real force_error=std::abs(network.normal_force_sum-old_force)/std::max(old_force,1e-10*reaction);
            const Real energy_error=std::abs(growth.contact_energy-old_energy)/std::max(old_energy,1e-22*reaction);
            Real nodal_force_change=0;
            for(size_t j=0;j<count;++j)nodal_force_change+=(bulk_force[j]-expected[j]).squaredNorm();
            nodal_force_change=std::sqrt(nodal_force_change)/std::max(reaction,Real(1e-30));
            const bool pass=corner_error<=1e-12&&opening_error<=1e-12&&force_error<=.01&&energy_error<=.01&&
                network.force_residual<=1e-8&&network.torque_residual<=1e-8&&network.power_residual<=1e-8&&
                new_step>old_step&&cb.commits==old_commits;
            std::ofstream audit(out/"affine_crack_contact_handoff.json");
            audit<<std::setprecision(17)<<"{\"status\":\""<<(pass?"PASSED":"FAILED")<<"\",\"step\":"<<time.steps
                <<",\"strain\":"<<u/o.diameter<<",\"reaction_N\":"<<reaction
                <<",\"elapsed_s\":"<<elapsed<<",\"old_crack_gap_m\":"<<old_gap<<",\"new_crack_gap_m\":"<<growth.max_gap
                <<",\"maximum_patch_opening_change_m\":"<<opening_error<<",\"maximum_corner_change_m\":"<<corner_error
                <<",\"old_penetration_m\":"<<old_penetration<<",\"new_penetration_m\":"<<growth.max_penetration
                <<",\"old_closed_contacts\":"<<old_active<<",\"new_closed_contacts\":"<<network.active_contacts
                <<",\"old_closed_normal_force_N\":"<<old_force<<",\"new_closed_normal_force_N\":"<<network.normal_force_sum
                <<",\"closed_force_relative_error\":"<<force_error<<",\"nodal_force_change_over_reaction\":"<<nodal_force_change
                <<",\"old_crack_contact_energy_J\":"<<old_energy<<",\"new_crack_contact_energy_J\":"<<growth.contact_energy
                <<",\"crack_contact_energy_relative_error\":"<<energy_error
                <<",\"force_balance_residual\":"<<network.force_residual<<",\"torque_balance_residual\":"<<network.torque_residual
                <<",\"contact_power_residual\":"<<network.power_residual
                <<",\"old_max_weight_l1\":"<<*std::max_element(old_weights.begin(),old_weights.end())<<",\"new_max_weight_l1\":"<<new_weight_max
                <<",\"old_crack_dt_s\":"<<old_dt<<",\"new_crack_dt_s\":"<<growth.contact_dt
                <<",\"old_combined_dt_s\":"<<old_step<<",\"new_combined_dt_s\":"<<new_step
                <<",\"material_steps_advanced\":"<<cb.commits-old_commits<<"}\n";
            if(!pass)throw std::runtime_error("Affine crack-contact handoff failed geometry, contact, or time-step verification");
            for(size_t j=0;j<count;++j)time.acceleration[j]=bulk_force[j]/mass[j];
        }
        if(affine_inner_interface&&recovery::checkpoint.read_version<10){
            const Real u=compression(stage_time),upper=contact.upper,lower=contact.lower,reaction=.5*(upper+lower);
            const Real old_dt=tie->stable_dt,old_energy=tie->energy,old_gap=tie->max_gap,crack_gap=growth.max_gap;
            const Real old_weight=si->maximum_weight_l1;
            const Real plate_dt=local_platen_step_bound?boundRows(nullptr,u):contact.stable_dt;
            const Real shell_dt=sb.finiteStrainTimeStep();
            auto stableStep=[&](Real interface_dt){
                Real rate=1/(core_dt*core_dt)+1/(shell_dt*shell_dt)+1/(sb.stabilizer_dt*sb.stabilizer_dt)
                    +1/(interface_dt*interface_dt)+1/(plate_dt*plate_dt);
                if(std::isfinite(growth.contact_dt))rate+=1/(growth.contact_dt*growth.contact_dt);
                return std::min({o.max_dt,1/std::sqrt(rate),cb.acoustic_dt,radial_dt});
            };
            const Field old_positions=tie->x_a,old_forces=tie->force_a;
            const size_t old_commits=cb.commits;
            inner_affine_positions.assign(sx,sx+sn);
            installAffineInner(old_positions,{});
            assembleShellAndConstraints(u);
            growth.statistics(growth.positions());
            Real position_change=0,force_change=0,force_scale=0,redistribution=0,new_residual=0;
            for(size_t i=0;i<old_positions.size();++i){
                position_change=std::max(position_change,(tie->x_a[i]-old_positions[i]).norm());
                force_change=std::max(force_change,(tie->force_a[i]-old_forces[i]).norm());
                force_scale=std::max(force_scale,old_forces[i].norm());
            }
            for(size_t j=0;j<count;++j){
                redistribution+=(bulk_force[j]-expected[j]).squaredNorm();
                new_residual+=bulk_force[j].squaredNorm();
            }
            const Real scale_N=std::max(reaction,Real(1e-30));
            redistribution=std::sqrt(redistribution)/scale_N;
            new_residual=std::sqrt(new_residual)/scale_N;
            const Real force_error=force_change/std::max(force_scale,1e-10*scale_N);
            const Real energy_error=std::abs(tie->energy-old_energy)/std::max(old_energy,1e-22*scale_N);
            const Real reaction_error=std::max(std::abs(contact.upper-upper),std::abs(contact.lower-lower))/scale_N;
            const Real old_step=stableStep(old_dt),new_step=stableStep(tie->stable_dt);
            const bool pass=position_change<=1e-12&&std::abs(growth.max_gap-crack_gap)<=1e-12&&force_error<=.01&&
                energy_error<=.01&&reaction_error<=.01&&tie->force_residual<=1e-8&&tie->torque_residual<=1e-8&&
                tie->power_residual<=1e-8&&new_step>old_step&&std::isfinite(new_residual)&&cb.commits==old_commits;
            // Reuse the usual fixed-displacement equilibration when the force
            // map change leaves a transient. The saved ramp clock stays fixed.
            inner_hold=pass&&new_residual>o.force_tolerance;
            if(inner_hold)checks=0;
            std::ofstream audit(out/"affine_inner_interface_handoff.json");
            audit<<std::setprecision(17)<<"{\"status\":\""<<(pass?"PASSED":"FAILED")<<"\",\"step\":"<<time.steps
                <<",\"strain\":"<<u/o.diameter<<",\"elapsed_s\":"<<elapsed<<",\"old_reaction_N\":"<<reaction
                <<",\"new_reaction_N\":"<<.5*(contact.upper+contact.lower)<<",\"reaction_relative_error\":"<<reaction_error
                <<",\"old_interface_gap_m\":"<<old_gap<<",\"new_interface_gap_m\":"<<tie->max_gap
                <<",\"maximum_interface_position_change_m\":"<<position_change
                <<",\"maximum_interface_force_change_N\":"<<force_change<<",\"interface_force_relative_error\":"<<force_error
                <<",\"old_interface_energy_J\":"<<old_energy<<",\"new_interface_energy_J\":"<<tie->energy
                <<",\"interface_energy_relative_error\":"<<energy_error
                <<",\"interface_force_balance_residual\":"<<tie->force_residual
                <<",\"interface_torque_balance_residual\":"<<tie->torque_residual
                <<",\"interface_power_residual\":"<<tie->power_residual
                <<",\"nodal_force_redistribution_over_reaction\":"<<redistribution<<",\"new_force_residual\":"<<new_residual
                <<",\"old_crack_gap_m\":"<<crack_gap<<",\"new_crack_gap_m\":"<<growth.max_gap
                <<",\"old_max_inner_weight_l1\":"<<old_weight<<",\"new_max_inner_weight_l1\":"<<si->maximum_weight_l1
                <<",\"inner_horizon_fallback_rows\":"<<si->extended_position_fallback_points
                <<",\"old_interface_dt_s\":"<<old_dt<<",\"new_interface_dt_s\":"<<tie->stable_dt
                <<",\"old_combined_dt_s\":"<<old_step<<",\"new_combined_dt_s\":"<<new_step
                <<",\"fixed_displacement_hold_required\":"<<(inner_hold?"true":"false")
                <<",\"material_steps_advanced\":"<<cb.commits-old_commits<<"}\n";
            if(!pass)throw std::runtime_error("Affine core-shell handoff failed force, work, gap or time-step verification");
            if(inner_hold)std::cout<<"AFFINE_INNER_REEQUILIBRATION_STARTED strain="<<u/o.diameter<<" residual="<<new_residual<<std::endl;
            for(size_t j=0;j<count;++j)time.acceleration[j]=bulk_force[j]/mass[j];
        }
        cohesive_force=growth.cohesiveForces(growth.positions());snapshot("resumed",compression(stage_time));
        std::ofstream receipt(out/"resume_validation.json");receipt<<std::setprecision(17)<<"{\"status\":\"PASSED\",\"step\":"<<time.steps<<",\"strain\":"<<compression(stage_time)/o.diameter<<",\"force_relative_error\":"<<error/std::max(scale,Real(1e-20))<<",\"core_commits\":"<<cb.commits<<",\"kinetic_energy_J\":"<<time.kinetic()<<",\"loading_history_reset\":false}\n";}
    save();
    if(column_handoff_only||affine_contact_handoff_only||inner_interface_handoff_only)return 0;

    while(stage<=o.load_steps){size_t prior_stage=stage,prior_cracks=growth.insertion_events;
        if(time.steps%50==0){for(size_t event=0;event<bank.size();++event){size_t k=growth.select();if(k==bank.size())break;
            surface_transfer_prepared=false;
            Real ratio=growth.strengthRatio(k);size_t prior_core_commits=cb.commits;bool inserted=growth.insert(k,false,refresh,eventEnergy,prepareTransfer,rollbackTransfer);if(cb.commits!=prior_core_commits)throw std::runtime_error("Topology event advanced frozen core material history");
            if(!inserted&&!growth.last_candidate_unresolved){snapshot("stopped_before_rejected_crack",compression(stage_time));save();throw std::runtime_error("Capsule physical insertion rejected: "+growth.last_rejection);}
            if(inserted){onset=true;cohesive_force=growth.cohesiveForces(growth.positions());events<<(fracture_contact_only?"strength_triggered_bond_deletion":"new_cohesive_crack")<<','<<time.steps<<','<<stage<<','<<compression(stage_time)/o.diameter<<','<<bank[k].id<<','<<ratio<<','<<sr.physical_broken_bonds<<','<<growth.max_gap<<'\n';events.flush();
                surface_audit<<time.steps<<','<<stage<<','<<compression(stage_time)/o.diameter<<','<<bank[k].id<<','<<(si?si->reassignedPoints():0)<<','<<so.reassignedPoints()<<','<<std::max(si?si->affine_error:Real(0),so.affine_error)<<','<<std::max(si?si->polynomial_error:Real(0),so.polynomial_error)<<','<<std::max(si?si->maximum_weight_l1:Real(0),so.maximum_weight_l1)<<','<<growth.last_insertion.fraction<<','<<core_dt<<','<<sb.finiteStrainTimeStep()<<','<<sb.stabilizer_dt<<','<<tie->stable_dt<<','<<contact.stable_dt<<','<<growth.contact_dt<<','<<so.position_fallback_points<<','<<so.extended_position_fallback_points<<'\n';surface_audit.flush();}
            else{events<<"unresolved_geometric_facet_no_material_cut"<<','<<time.steps<<','<<stage<<','<<compression(stage_time)/o.diameter<<','<<bank[k].id<<','<<ratio<<','<<sr.physical_broken_bonds<<','<<growth.max_gap<<'\n';events.flush();}
        }}
        // A cascade changes topology at one fixed physical time. Preserve
        // every event in CSV and publish its final state once, avoiding many
        // redundant full-core/shell files at the same timestamp.
        if(growth.insertion_events!=prior_cracks)snapshot("crack_"+std::to_string(growth.insertion_events),compression(stage_time));
        const Real shell_dt=sb.finiteStrainTimeStep();
        auto allowedStep=[&](Real platen_dt){
            Real rate=1/(core_dt*core_dt)+1/(shell_dt*shell_dt)+1/(sb.stabilizer_dt*sb.stabilizer_dt)
                +1/(tie->stable_dt*tie->stable_dt)+1/(platen_dt*platen_dt);
            if(std::isfinite(growth.contact_dt))rate+=1/(growth.contact_dt*growth.contact_dt);
            Real bounded=std::min({o.max_dt,1/std::sqrt(rate),cb.acoustic_dt,radial_dt});
            return (!column_hold&&!inner_hold&&stage_time<o.ramp_duration)?std::min(bounded,o.ramp_duration-stage_time):bounded;
        };
        const Real all_rows_dt=allowedStep(contact.stable_dt),old_u=compression(stage_time);
        Real dt=local_platen_step_bound?allowedStep(boundRows(nullptr,old_u)):all_rows_dt;

        Field old_d=time.displacement,old_internal=internal,old_shell,vhalf,qfree,next_d;
        if(!fracture_contact_only){old_shell.assign(sx,sx+sn);qfree.resize(sn);}
        const Real oldReaction=.5*(contact.upper+contact.lower),oldK=time.kinetic(),psiold=growth.stored;
        Real dampedK=0,factor=1;
        CohesiveSurfaceEvolution::Result cohesive;
        for(size_t retry=0;;++retry){
            if(!std::isfinite(dt)||dt<1e-15||time.steps>=o.max_steps)
                throw std::runtime_error("Coupled physical time/count bound exhausted");
            factor=std::exp(-.5*o.damping_rate*dt);
            vhalf=time.velocity;dampedK=0;
            for(size_t j=0;j<count;++j){vhalf[j]*=factor;dampedK+=.5*mass[j]*vhalf[j].squaredNorm();vhalf[j]+=.5*dt*bulk_force[j]/mass[j];}
            next_d=old_d;
            for(size_t j=0;j<count;++j)next_d[j]+=dt*vhalf[j];
            if(!fracture_contact_only){
                for(size_t j=0;j<sn;++j)qfree[j]=sX[j]+next_d[cn+j];
                try{
                    cohesive=growth.trialDrift(qfree,old_shell,inverse_shell_mass,dt);
                }catch(const std::runtime_error&){
                    if(dt<=all_rows_dt)throw;
                    dt=all_rows_dt;continue;
                }
            }
            if(dt<=all_rows_dt)break;
            // The drift follows a chord between the two trace positions,
            // and the planes move inwards. Include all rows that can touch
            // either plane in this step. Reject before committing histories.
            if(!fracture_contact_only)for(size_t j=0;j<sn;++j)next_d[cn+j]=cohesive.positions[j]-sX[j];
            const Real checked_dt=allowedStep(boundRows(&next_d,compression(stage_time+dt)));
            if(dt<=checked_dt)break;
            dt=retry>=2?std::min(all_rows_dt,checked_dt):checked_dt;
        }
        Wd+=oldK-dampedK;
        if(!fracture_contact_only){
            growth.acceptDrift(cohesive);cohesive_force=cohesive.force;
            Real velocity_defect=0,surface_work=0;for(size_t j=0;j<sn;++j){next_d[cn+j]=cohesive.positions[j]-sX[j];Vecd after=vhalf[cn+j]+dt*cohesive.force[j]/sm[j];velocity_defect+=.5*sm[j]*(after-vhalf[cn+j]).squaredNorm();surface_work-=cohesive.force[j].dot(cohesive.positions[j]-old_shell[j]);vhalf[cn+j]=after;}
            Walgorithmic+=velocity_defect+surface_work-(cohesive.stored-psiold)-cohesive.dissipation_increment;
        }
        time.displacement=next_d;forceAt(time.displacement,vhalf,compression(stage_time+dt),true);
        Real undampedK=0;for(size_t j=0;j<count;++j){time.velocity[j]=vhalf[j]+.5*dt*bulk_force[j]/mass[j];undampedK+=.5*mass[j]*time.velocity[j].squaredNorm();time.velocity[j]*=factor;time.acceleration[j]=(bulk_force[j]+(j>=cn?cohesive_force[j-cn]:Vecd::Zero()))/mass[j];}
        Wd+=undampedK-time.kinetic();++time.steps;setPositions(time.displacement,time.velocity);Real u=compression(stage_time+dt),mean=.5*(contact.upper+contact.lower);Wp+=.5*(oldReaction+mean)*(u-old_u);for(size_t j=0;j<cn;++j)Wcore+=.5*(old_internal[j]+internal[j]).dot(time.displacement[j]-old_d[j]);if(inner_hold)inner_hold_elapsed+=dt;else stage_time+=dt;elapsed+=dt;
        if(!fracture_contact_only||time.steps%200==0||(!fully_failed&&onset))growth.statistics(growth.positions());
        if(!fully_failed&&growth.failed_quadrature){fully_failed=true;events<<(fracture_contact_only?"first_contact_only_crack":"first_full_cohesive_failure")<<','<<time.steps<<','<<stage<<','<<u/o.diameter<<",0,0,"<<sr.physical_broken_bonds<<','<<growth.max_gap<<'\n';events.flush();snapshot("first_full_failure",u);}
        if(time.steps%200==0){Real residual=0,minc=1,mins=1,newK=time.kinetic();for(size_t j=0;j<count;++j){Vecd f=bulk_force[j];if(j>=cn)f+=cohesive_force[j-cn];residual+=f.squaredNorm();}residual=std::sqrt(residual)/std::max(mean,Real(1e-20));for(const auto&F:cb.F)minc=std::min(minc,F.determinant());for(size_t j=0;j<sn;++j)mins=std::min(mins,sF[j].determinant());Real imbalance=std::abs(contact.upper-contact.lower)/std::max(mean,Real(1e-20)),kr=newK/std::max(mean*u,Real(1e-30));Real budget=Wp+interface_switch_work-Wcore-(recoverable()-U0)-(growth.dissipated-growth.initial_dissipation)-(newK-K0)-Wd-Walgorithmic-growth.stabilizer_transfer_loss+growth.insertion_defect;
            auto row=[&](std::ostream&f){f<<time.steps<<','<<stage<<','<<elapsed<<','<<u/o.diameter<<','<<u<<','<<dt<<','<<contact.upper<<','<<contact.lower<<','<<residual<<','<<imbalance<<','<<kr<<','<<minc<<','<<mins<<','<<sr.physical_broken_bonds<<','<<growth.patches.size()<<','<<growth.failed_quadrature<<','<<growth.max_gap<<','<<growth.stored<<','<<growth.dissipated-growth.initial_dissipation<<','<<sb.energy()<<','<<tie->energy<<','<<contact.energy<<','<<growth.contact_energy<<','<<std::max(contact.max_penetration,growth.max_penetration)<<','<<Wp<<','<<Wcore<<','<<Wd<<','<<Walgorithmic<<','<<growth.stabilizer_transfer_loss<<','<<growth.insertion_defect<<','<<budget<<'\n';f.flush();};row(progress);if(time.steps%1000==0)std::cout<<"PROGRESS step="<<time.steps<<" stage="<<stage<<" strain="<<u/o.diameter<<" force="<<mean<<" residual="<<residual<<" cracks="<<growth.patches.size()<<" gap="<<growth.max_gap<<" dt="<<dt<<std::endl;
            bool good=(column_hold||inner_hold||stage_time>=o.ramp_duration)&&mean>1e-10&&residual<=o.force_tolerance&&imbalance<=o.balance_tolerance&&kr<=o.kinetic_tolerance;
            checks=good?checks+1:0;
            if(checks>=3){
                if(column_hold||inner_hold){
                    auto file=out/(inner_hold?"affine_inner_reequilibration.json":"column_reequilibration.json");
                    std::ofstream report(file);report<<std::setprecision(17)
                        <<"{\"status\":\"PASSED\",\"step\":"<<time.steps<<",\"strain\":"<<u/o.diameter
                        <<",\"hold_elapsed_s\":"<<(inner_hold?inner_hold_elapsed:stage_time)<<",\"reaction_N\":"<<mean<<",\"force_residual\":"<<residual
                        <<",\"kinetic_ratio\":"<<kr<<",\"dt_s\":"<<dt<<",\"crack_facets\":"<<growth.patches.size()<<"}\n";
                    snapshot(inner_hold?"inner_rebalanced":"column_rebalanced",u);
                    const char*event=inner_hold?"AFFINE_INNER_REEQUILIBRATION_ACCEPTED":"FIXED_COLUMN_REEQUILIBRATION_ACCEPTED";
                    if(column_hold)stage_time=0;
                    column_hold=inner_hold=false;inner_hold_elapsed=0;checks=0;
                    std::cout<<event<<" force="<<mean<<" dt="<<dt<<std::endl;save();
                }else{
                    row(eq);snapshot(std::to_string(stage),u);
                    std::cout<<"ACCEPTED stage="<<stage<<" force="<<mean<<" cracks="<<growth.patches.size()<<std::endl;
                    ++stage;stage_time=0;checks=0;
                }
            }
        }
        if(inner_hold&&inner_hold_elapsed>o.max_stage_duration){save();throw std::runtime_error("Affine inner-interface hold did not equilibrate");}
        if(stage_time>o.max_stage_duration)throw std::runtime_error("Physical capsule stage did not equilibrate");if(time.steps%200==0&&recovery::checkpoint.due(time.steps,stage!=prior_stage))save();if(recovery::checkpoint.stop_step&&time.steps>=recovery::checkpoint.stop_step){save();return 86;}if(recovery::checkpoint.stop_after_equilibrium&&stage!=prior_stage){save();return 86;}
    }
    save();std::ofstream summary(out/"summary.json");summary<<"{\"status\":\"COMPLETE_PHYSICAL_FRACTURE_MECHANICAL_HISTORY\",\"accepted_points\":"<<o.load_steps<<",\"steps\":"<<time.steps<<",\"automatic_crack_facets\":"<<growth.patches.size()<<",\"broken_bond_pairs\":"<<sr.physical_broken_bonds<<",\"fracture_contact_only\":"<<(fracture_contact_only?"true":"false")<<",\"maximum_current_crack_gap_m\":"<<growth.max_gap<<",\"full_cohesive_failure\":"<<(!fracture_contact_only&&fully_failed?"true":"false")<<",\"parameters_changed\":false}\n";return 0;
}

int main(int argc,char**argv)
{
    try{recovery::checkpoint.schema_version=4;std::string cells;PeriExecution::LegacyGpuSelection legacy{};std::vector<char*>args{argv[0]};for(int i=1;i<argc;++i){std::string a=argv[i];if(a.rfind("--snsni-cells=",0)==0)cells=a.substr(14);else if(a.rfind("--fracture-facets=",0)==0)fracture_bank_input=a.substr(18);else if(a=="--column-tied-interface")use_column_interface=true;else if(a=="--column-interface-handoff-only")column_handoff_only=true;else if(a=="--local-platen-step-bound")local_platen_step_bound=true;else if(a=="--fracture-contact-only")fracture_contact_only=true;else if(a=="--affine-crack-contact")affine_crack_contact=true;else if(a=="--affine-contact-handoff-only")affine_contact_handoff_only=true;else if(a=="--affine-inner-interface")affine_inner_interface=true;else if(a=="--inner-interface-handoff-only")inner_interface_handoff_only=true;else if(a=="--cuda-crack-projection")legacy.crack_projection=true;else if(a=="--cuda-bonded-interface")legacy.bonded_interface=true;else if(recovery::checkpoint.take(a)){}else args.push_back(argv[i]);}
        if(use_column_interface){recovery::checkpoint.schema_version=5;recovery::checkpoint.compatible_read_version=4;}
        if(fracture_contact_only){
            recovery::checkpoint.compatible_read_version=recovery::checkpoint.schema_version;
            recovery::checkpoint.schema_version+=2;
        }
        if(affine_crack_contact){
            if(!fracture_contact_only)throw std::invalid_argument("Affine crack contact requires contact-only fracture");
            recovery::checkpoint.compatible_read_version=recovery::checkpoint.schema_version;
            recovery::checkpoint.schema_version+=2;
        }
        if(affine_inner_interface){
            if(!affine_crack_contact||use_column_interface)
                throw std::invalid_argument("Affine inner interface requires the affine crack-contact model with dynamic shell-side traces");
            recovery::checkpoint.compatible_read_version=recovery::checkpoint.schema_version;
            recovery::checkpoint.schema_version=10;
        }
        if(inner_interface_handoff_only&&(!affine_inner_interface||recovery::checkpoint.resume.empty()))
            throw std::invalid_argument("Inner interface handoff audit requires a saved state and --affine-inner-interface");
        if(affine_contact_handoff_only&&(!affine_crack_contact||recovery::checkpoint.resume.empty()))
            throw std::invalid_argument("Affine contact handoff audit requires a saved state and --affine-crack-contact");
        if(column_handoff_only&&(!use_column_interface||recovery::checkpoint.resume.empty()))
            throw std::invalid_argument("Column handoff audit requires the fixed-column interface and an accepted checkpoint");
        if(cells.empty())throw std::invalid_argument("--snsni-cells is required");HybridOptions options;options.parseHybrid(int(args.size()),args.data());
        if(local_platen_step_bound&&options.plate_friction!=0)
            throw std::invalid_argument("Local platen step bound requires frictionless analytic planes");
        const PeriExecution execution(options.execution,legacy);
        tbb::task_arena arena(int(options.threads));return arena.execute([&]{return runSnsniCapsule(options,cells,execution);});
    }catch(const std::exception&e){std::cerr<<"COUPLED_ERROR "<<e.what()<<std::endl;return 1;}
}
