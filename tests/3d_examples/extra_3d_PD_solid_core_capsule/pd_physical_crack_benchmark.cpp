#include "pd_snsni_capsule_fixture.h"
#include "pd_broken_bond_energy.h"
#include "pd_fracture_visible_surface.h"
#include "pd_crack_face_contact.h"
#include "pd_cohesive_crack_law.h"
#include "pd_fracture_topology_checkpoint.h"
#include <Eigen/Geometry>
#include <filesystem>
#include <iomanip>
using namespace SPH;using namespace SPH::solid_dynamics;using namespace SPH::peridynamics;
namespace SPH {
class PhysicalCrackSlab;
template<>class ParticleGenerator<BaseParticles,PhysicalCrackSlab>:public ParticleGenerator<BaseParticles>{
    Real dp,t;int nx,ny;
  public:ParticleGenerator(SPHBody&b,BaseParticles&p,Real d,Real h,int n,int m):ParticleGenerator<BaseParticles>(b,p),dp(d),t(h),nx(n),ny(m){}
    void prepareGeometricData()override{ThicknessQuadrature rule(3);for(size_t k=0;k<3;++k)for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)addPositionAndVolumetricMeasure(Vecd((i+.5-nx/2.)*dp,(j+.5-ny/2.)*dp,.5*t*rule.xi[k]),dp*dp*.5*t*rule.weight[k]);}
};}
void writeCloud(const std::filesystem::path&p,const StdVec<Vecd>&x,const StdVec<Real>&d){
 std::ofstream f(p);f<<std::setprecision(17)<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\"><PolyData><Piece NumberOfPoints=\""<<x.size()<<"\" NumberOfVerts=\""<<x.size()<<"\"><Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";for(const auto&a:x)f<<a.transpose()<<' ';f<<"</DataArray></Points><PointData><DataArray type=\"Float64\" Name=\"BrokenBondFraction\" format=\"ascii\">";for(Real a:d)f<<a<<' ';f<<"</DataArray></PointData><Verts><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">";for(size_t i=0;i<x.size();++i)f<<i<<' ';f<<"</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">";for(size_t i=0;i<x.size();++i)f<<i+1<<' ';f<<"</DataArray></Verts></Piece></PolyData></VTKFile>\n";
}
void writeFace(const std::filesystem::path&p,const StdVec<Vecd>&x){std::ofstream f(p);f<<std::setprecision(17)<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\"><PolyData><Piece NumberOfPoints=\"4\" NumberOfPolys=\"2\"><Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";for(const auto&a:x)f<<a.transpose()<<' ';f<<"</DataArray></Points><Polys><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">0 1 2 0 2 3</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">3 6</DataArray></Polys></Piece></PolyData></VTKFile>\n";}
#include "pd_automatic_growth_benchmark.h"
int main(int argc,char**argv){bool published=false;try{
 if(argc!=2)throw std::invalid_argument("Output directory required");std::filesystem::path out(argv[1]);if(std::filesystem::exists(out)&&!std::filesystem::is_empty(out))throw std::runtime_error("Benchmark output must be new or empty");std::filesystem::create_directories(out/"frames");
 {std::ofstream pending(out/"acceptance.json");pending<<"{\"geometry_status\":\"RUNNING\",\"capsule_restart_allowed\":false}\n";if(!pending)throw std::runtime_error("Cannot publish precrack acceptance start");published=true;}
 const Real E=2.151e9,nu=.3,ft=46.4e6,Gc=100,dp=25e-6,t=5e-6,K=2.8491556e15;const int nx=16,ny=16;
 SPHSystem system(BoundingBoxd(Vecd::Constant(-.002),Vecd::Constant(.002)),dp,8);PDBody body(system,makeShared<DefaultShape>("RealPrecrackedNOSBLayer"));body.defineAdaptation<PDAdaptation>(1.5075,1);body.defineMatterMaterial<HughesWingetSolid>(1200,E,nu);body.generateParticles<BaseParticles,PhysicalCrackSlab>(dp,t,nx,ny);
 InnerRelation inner(body);SimpleDynamics<PDTimeStepInitialization> init(body,makeShared<Gravity>(Vecd::Zero()));SimpleDynamics<NosbPDFirstStep> first(body);InteractionWithUpdate<NosbPDShapeMatrix> shape(inner);auto&p=body.getBaseParticles();size_t count=p.TotalRealParticles(),nc=count/3;
 StdVec<StdVec<size_t>>columns(nc);StdVec<Vecd>directors(nc,Vecd::UnitZ());for(size_t c=0;c<nc;++c)for(size_t k=0;k<3;++k)columns[c].push_back(k*nc+c);ThicknessQuadrature rule(3);LayeredParticleGeometry geometry(p,columns,directors,rule.xi);PDLayeredFractureReference ref(inner,geometry);
 InteractionWithUpdate<NosbPDQuadraticSecondStep> second(inner,ref);InteractionDynamics<NosbPDQuadraticThirdStep> third(inner,ref);InteractionDynamics<ProjectedPairStabilization> unused_hg(inner,.5,.1);InteractionDynamics<PairNumericalDampingforPD> unused_damping(inner);
 system.initializeSystemCellLinkedLists();system.initializeSystemConfigurations();shape.exec();ref.initialize();auto intact_topology=ref.captureTopology();size_t broken=ref.cutPlane(Vecd::Zero(),Vecd::UnitX());if(!broken)throw std::runtime_error("Precrack cut no physical bonds");
 PDBrokenBondEnergy bulk(inner,ref,E,nu,.5,body.getSPHAdaptation().getKernel()->CutOffRadius());bulk.rebuild();State s(body);SolidForcePort port(body,SolidForcePort::Storage::PDAccelerationPrior);
 auto*acc=p.getVariableDataByName<Vecd>("PDAcceleration");auto*bd=p.getVariableDataByName<Real>("BrokenBondFraction");StdVec<Vecd>master_points={Vecd(0,-ny*dp/2,-t/2),Vecd(0,ny*dp/2,-t/2),Vecd(0,ny*dp/2,t/2),Vecd(0,-ny*dp/2,t/2)},slave_points;StdVec<Real>areas;
 for(int j=0;j<8;++j)for(size_t k=0;k<3;++k){slave_points.push_back(Vecd(0,(j+.5-4)*dp,.5*t*rule.xi[k]));areas.push_back(dp*.5*t*rule.weight[k]);}
 auto anchors=[&](const StdVec<Vecd>&points,int side){StdVec<size_t>result;for(const auto&q:points){size_t a=count;Real best=std::numeric_limits<Real>::infinity();for(size_t i=0;i<count;++i)if((s.x0[i][0]>0?1:-1)==side){Real d=(s.x0[i]-q).squaredNorm();if(d<best){best=d;a=i;}}if(a==count)throw std::runtime_error("Missing crack-face owner");result.push_back(a);}return result;};
 VisibleCrackSurface master(port,master_points,anchors(master_points,-1),ref,dp,3,2),slave(port,slave_points,anchors(slave_points,1),ref,dp,3,2),plus_corners(port,master_points,anchors(master_points,1),ref,dp,3,2);
 BilinearCrackSurface minus_geometry(master,slave_points),plus_geometry(plus_corners,slave_points);
 CrackSurfaceContactNetwork contact({{&minus_geometry,&areas,1,0},{&plus_geometry,&areas,-1,0}},K,dp);

 Real volume=0,area=0;for(size_t i=0;i<count;++i)volume+=s.measure[i];for(Real a:areas)area+=a;
 auto evaluate=[&](){for(size_t i=0;i<count;++i)bulk.evaluate(i);port.clear();for(size_t i=0;i<count;++i){acc[i].setZero();bulk.forces(i);}master.ensureCurrent();slave.ensureCurrent();contact.apply();for(size_t i=0;i<count;++i)acc[i]+=port.external[i];return bulk.energy()+contact.energy;};
 auto set=[&](Real gap,Real slide,const Matd&Q){for(size_t i=0;i<count;++i){Vecd jump=Vecd((s.x0[i][0]>0?.5:-.5)*gap,s.x0[i][0]>0?slide:0,0);s.x[i]=Q*(s.x0[i]+jump);s.v[i].setZero();}};
 StdVec<int>components(count,-1);int component_count=0;for(size_t start=0;start<count;++start)if(components[start]<0){StdVec<size_t>queue{start};components[start]=component_count;for(size_t k=0;k<queue.size();++k){size_t i=queue[k];for(size_t n=0;n<ref.rows[i].size();++n)if(ref.live(i,n)){size_t j=ref.rows[i][n].j;if(components[j]<0){components[j]=component_count;queue.push_back(j);}}}++component_count;}
 size_t bad_bonds=0;for(size_t i=0;i<count;++i)for(size_t n=0;n<ref.rows[i].size();++n){const auto&e=ref.rows[i][n];bad_bonds+=ref.live(i,n)&&(s.x0[i][0]*s.x0[e.j][0]<0);if(!ref.live(i,n)&&e.gradient.norm()!=0)throw std::runtime_error("A broken bond retains a material gradient");}
 // The original point-level Gc history may be retained through a fully
 // decohered continuum-to-discrete transition, but only if a separate growth
 // controller verifies the actual transition energy/traction. Here check the
 // live residual itself and ensure rebuilding does not reset fracture history.
 PDEnergyCrackBand live_damage(inner,ref,E,nu,ft,Gc,.5,body.getSPHAdaptation().getKernel()->CutOffRadius(),true,true);live_damage.initialize();
 for(size_t i=0;i<count;++i){live_damage.history[i]={live_damage.fracture.kf,1,Gc/live_damage.fracture.width};live_damage.crack_normal[i]=Vecd::UnitX();}
 set(.0006,0,Matd::Identity());live_damage.initialize();bool damage_history_preserved=true;Real live_damage_force=0;
 for(size_t i=0;i<count;++i){damage_history_preserved=damage_history_preserved&&live_damage.history[i].damage==1&&live_damage.history[i].kappa==live_damage.fracture.kf&&live_damage.history[i].dissipation==Gc/live_damage.fracture.width;live_damage.evaluate(i,false);}
 for(size_t i=0;i<count;++i){acc[i].setZero();live_damage.forces(i);live_damage_force=std::max(live_damage_force,(s.mass[i]*acc[i]).norm()/(E*area));}
 Real live_damage_energy=std::abs(live_damage.statistics().stored)/(E*volume);
 // Cache safety belongs to this same specimen: deliberately miss several
 // topology epochs and verify no force/contact path can use an old operator.
 auto must_reject=[](auto fn){try{fn();}catch(const std::exception&){return true;}return false;};
 auto saved_topology=ref.captureTopology();bool severed=false;
 for(size_t i=0;i<count&&!severed;++i)for(size_t n=0;n<ref.rows[i].size()&&!severed;++n)
  if(ref.live(i,n)&&geometry.column_of[i]!=geometry.column_of[ref.rows[i][n].j])severed=ref.severBondPair(i,n);
 if(!severed)throw std::runtime_error("No live tangent bond for topology guard verification");
 bool pending_guard=must_reject([&](){bulk.forces(0);})&&must_reject([&](){live_damage.forces(0);})&&must_reject([&](){contact.apply();});
 ref.restoreTopology(saved_topology);ref.restoreTopology(saved_topology);
 bool stale_guard=must_reject([&](){bulk.evaluate(0);})&&must_reject([&](){live_damage.evaluate(0,false);})&&must_reject([&](){master.position(0);})&&must_reject([&](){contact.apply();});
 master.rebuildChanged();slave.rebuildChanged();plus_corners.rebuildChanged();minus_geometry.rebuild();plus_geometry.rebuild();bulk.rebuild();live_damage.initialize();
 bool topology_guards=pending_guard&&stale_guard&&master.isCurrent()&&slave.isCurrent();
 if(!topology_guards)throw std::runtime_error("Stale/pending physical fracture cache was accepted");
 // This is a topology/history codec check on the SAME specimen. It does not
 // turn this test into a complete mechanical capsule restart.
 using Archive=FractureTopologyCheckpoint;std::vector<uint64_t>point_ids(count);
 for(size_t i=0;i<count;++i)point_ids[i]=i;
 const std::string archive_tag="precrack-codec-E2151MPa-nu0.3-ft46.4MPa-Gc100-K28491556e8-v1";
 Archive::Histories histories;Archive::IMEXBlock imex_block;imex_block.name="codec_only_initially_rigid_law";imex_block.ft=ft;imex_block.Gc=Gc;imex_block.history_ids={101};imex_block.history={{0}};histories.imex.push_back(imex_block);
 Archive::FiniteKBlock finite_block;finite_block.name="codec_only_finite_K_law";finite_block.ft=ft;finite_block.Gc=Gc;finite_block.K=K;finite_block.history_ids={201};BilinearCohesiveCrackLaw codec_law(ft,Gc,K);finite_block.history={codec_law.precrackedHistory()};histories.finite_k.push_back(finite_block);
 auto archive=Archive::encode(ref,archive_tag,point_ids,histories),corrupted=archive;corrupted.at(corrupted.size()/2)^=1;
 size_t before_bad_load=ref.topologyVersion();
 bool checksum_rejected=must_reject([&](){Archive::load(ref,corrupted,archive_tag,point_ids,histories);});
 bool identity_rejected=must_reject([&](){Archive::load(ref,archive,archive_tag+"-wrong",point_ids,histories);});
 auto wrong_layout=histories;wrong_layout.finite_k[0].name+="-wrong";
 bool layout_rejected=must_reject([&](){Archive::load(ref,archive,archive_tag,point_ids,wrong_layout);});
 bool rejected_without_mutation=ref.topologyVersion()==before_bad_load;
 Archive::load(ref,archive,archive_tag,point_ids,histories);
 bool checkpoint_roundtrip=archive==Archive::encode(ref,archive_tag,point_ids,histories)&&checksum_rejected&&identity_rejected&&layout_rejected&&rejected_without_mutation;
 if(!checkpoint_roundtrip)throw std::runtime_error("Physical fracture checkpoint roundtrip/rejection failed");
 std::ofstream archived(out/"fracture_topology_checkpoint.bin",std::ios::binary);archived.write(reinterpret_cast<const char*>(archive.data()),archive.size());if(!archived)throw std::runtime_error("Cannot persist fracture topology benchmark archive");archived.close();
 master.rebuildChanged();slave.rebuildChanged();plus_corners.rebuildChanged();minus_geometry.rebuild();plus_geometry.rebuild();bulk.rebuild();live_damage.initialize();


 struct Mode{const char*name;Real gap,slide;Real angle;};StdVec<Mode>modes={{"closed_initial",0,0,0},{"opening",20e-6,0,0},{"wide_open",.0006,0,0},{"open_sliding",20e-6,2*dp,0},{"closure",-2e-9,0,0},{"closed_sliding",-2e-9,2*dp,0},{"rotated_closed_sliding",-2e-9,2*dp,.37},{"reopening",20e-6,2*dp,.37}};
 std::ofstream csv(out/"opening_closure_sliding.csv");csv<<std::setprecision(17)<<"frame,mode,gap_m,slip_m,bulk_J,contact_J,normal_N,expected_normal_N,dt_bound_s\n";std::ofstream pvd(out/"precrack.pvd");pvd<<"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"><Collection>\n";
 Real max_bulk=0,max_contact_error=0,max_fd=0,max_balance=0,max_torque=0,max_slip=0,max_gap_error=0;size_t outside=0;
 for(size_t frame=0;frame<modes.size();++frame){auto m=modes[frame];Matd Q=Eigen::AngleAxis<Real>(m.angle,Vecd(1,2,3).normalized()).toRotationMatrix();set(m.gap,m.slide,Q);evaluate();max_bulk=std::max(max_bulk,std::abs(bulk.energy())/(E*volume));Real expected=K*area*std::max(-m.gap,Real(0));max_contact_error=std::max(max_contact_error,std::abs(contact.normal_force_sum-expected)/std::max(expected,E*area));max_slip=std::max(max_slip,contact.tangential_force_sum/(E*area));max_balance=std::max(max_balance,contact.force_residual);max_torque=std::max(max_torque,contact.torque_residual);outside+=contact.outside_projections;
  Real mean_gap=0;Vecd normal=Q*Vecd::UnitX();for(size_t q=0;q<slave.rows.size();++q)mean_gap+=normal.dot(slave.position(q)-master.position(0))/slave.rows.size();max_gap_error=std::max(max_gap_error,std::abs(mean_gap-m.gap)/dp);
  for(size_t i=0;i<count;++i)s.v[i]=Vecd(.012*std::sin(5*s.x0[i][0]/dp),.008*std::cos(3*s.x0[i][1]/dp),.004*std::sin(2*s.x0[i][2]/t));evaluate();Real power=0;for(size_t i=0;i<count;++i)power+=s.mass[i]*acc[i].dot(s.v[i]);
  if(m.gap<0){Real eps=1e-10;for(size_t i=0;i<count;++i)s.x[i]+=eps*s.v[i];Real ep=evaluate();for(size_t i=0;i<count;++i)s.x[i]-=2*eps*s.v[i];Real em=evaluate();max_fd=std::max(max_fd,std::abs((ep-em)/(2*eps)+power)/std::max(std::abs(power),Real(1e-12)));}
  set(m.gap,m.slide,Q);evaluate();Real dt=std::min({bulk.finiteStrainTimeStep(),bulk.stabilizer_dt,contact.stable_dt});csv<<frame<<','<<m.name<<','<<m.gap<<','<<m.slide<<','<<bulk.energy()<<','<<contact.energy<<','<<contact.normal_force_sum<<','<<expected<<','<<dt<<'\n';
  StdVec<Vecd>left,right;StdVec<Real>ld,rd;for(size_t i=0;i<count;++i){if(s.x0[i][0]<0){left.push_back(s.x[i]);ld.push_back(bd[i]);}else{right.push_back(s.x[i]);rd.push_back(bd[i]);}}
  auto name=std::to_string(frame);writeCloud(out/"frames"/("left_"+name+".vtp"),left,ld);writeCloud(out/"frames"/("right_"+name+".vtp"),right,rd);StdVec<Vecd>lm,rm;for(const auto&q:master_points){lm.push_back(Q*(q+Vecd(-.5*m.gap,0,0)));rm.push_back(Q*(q+Vecd(.5*m.gap,m.slide,0)));}writeFace(out/"frames"/("minus_face_"+name+".vtp"),lm);writeFace(out/"frames"/("plus_face_"+name+".vtp"),rm);
  std::ofstream vtm(out/"frames"/("frame_"+name+".vtm"));vtm<<"<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\"><vtkMultiBlockDataSet>";int block=0;for(auto part:{"left","right","minus_face","plus_face"})vtm<<"<DataSet index=\""<<block++<<"\" name=\""<<part<<"\" file=\""<<part<<'_'<<name<<".vtp\"/>";vtm<<"</vtkMultiBlockDataSet></VTKFile>\n";pvd<<"<DataSet timestep=\""<<frame<<"\" part=\"0\" file=\"frames/frame_"<<name<<".vtm\"/>\n";
 }
 pvd<<"</Collection></VTKFile>\n";
 // Same specimen, fixed displacement hold with physical masses and NO damping.
 set(-2e-9,dp,Matd::Identity());StdVec<int>fixed(count);for(size_t i=0;i<count;++i){size_t c=i%nc,ix=c/ny,iy=c%ny;fixed[i]=ix==0||ix==nx-1||iy==0||iy==ny-1;s.v[i].setZero();}
 Real e0=evaluate(),dt=.1*std::min({bulk.finiteStrainTimeStep(),bulk.stabilizer_dt,contact.stable_dt}),drift=0;for(int k=0;k<200;++k){for(size_t i=0;i<count;++i)if(!fixed[i]){s.v[i]+=.5*dt*acc[i];s.x[i]+=dt*s.v[i];}Real ep=evaluate(),kinetic=0;for(size_t i=0;i<count;++i)if(!fixed[i]){s.v[i]+=.5*dt*acc[i];kinetic+=.5*s.mass[i]*s.v[i].squaredNorm();}drift=std::max(drift,std::abs((ep+kinetic)/e0-1));if(dt>std::min({bulk.finiteStrainTimeStep(),bulk.stabilizer_dt,contact.stable_dt}))throw std::runtime_error("Explicit test exceeded current physical time bound");}
 // The SAME crack faces carry one initially cohesive ligament to verify the
 // independent ft/Gc law. This is NOT a test of dynamic topology insertion.
 BilinearCohesiveCrackLaw law(ft,Gc,K);BilinearCohesiveCrackLaw::History hist;Real work=0,last_t=0,last_gap=0,peak=0;const int steps=12000;
 for(int k=1;k<=steps;++k){Real gap=law.delta0*std::min(k,1000)/1000.+(k>1000?(law.delta_f-law.delta0)*(k-1000)/(steps-1000):0);auto r=law.evaluate(Vecd(gap,0,0),Vecd::UnitX(),hist,true);work+=.5*(last_t+r.traction[0])*(gap-last_gap);last_t=r.traction[0];last_gap=gap;peak=std::max(peak,r.traction[0]);}
 Real gc_error=std::abs(work/Gc-1),strength_error=std::abs(peak/ft-1);auto full=law.evaluate(Vecd(2*law.delta_f,dp,0),Vecd::UnitX(),hist,true);Real fully_failed_traction=full.traction.norm()/ft;
 bool geometry_pass=checkpoint_roundtrip&&topology_guards&&component_count==2&&bad_bonds==0&&broken>0&&max_bulk<1e-12&&max_contact_error<1e-9&&max_slip<1e-10&&max_gap_error<1e-9&&max_fd<1e-5&&max_balance<1e-10&&max_torque<1e-9&&outside==0&&drift<1e-3&&gc_error<1e-5&&strength_error<1e-10&&fully_failed_traction==0&&damage_history_preserved&&live_damage_energy<1e-12&&live_damage_force<1e-10;
 {std::ofstream stage(out/"precrack_stage_acceptance.json");stage<<std::setprecision(17)<<"{\"status\":\""<<(geometry_pass?"PASSED":"FAILED")<<"\",\"broken_bond_pairs\":"<<broken<<",\"force_energy_relative_error\":"<<max_fd<<",\"force_balance\":"<<max_balance<<",\"torque_balance\":"<<max_torque<<",\"explicit_energy_relative_drift\":"<<drift<<"}\n";}
 if(!geometry_pass)throw std::runtime_error("Precrack geometry/energy acceptance failed");
  // Reproduce the capsule's failed boundary-owner update on THIS same slab.
  // A reference boundary point can be across a new crack from its original
  // nearest-node search anchor. Fixed crack-side ownership must still reject it.
  {
   auto saved=ref.captureTopology();StdVec<Vecd>saved_x(s.x,s.x+count),saved_v(s.v,s.v+count);
   ref.restoreTopology(intact_topology);bulk.rebuild();
   for(size_t i=0;i<count;++i){const auto&x=s.x0[i];s.x[i]=x+Vecd(.002*x[0]+.002*x[0]*x[1]/(nx*dp),-.001*x[1],.001*x[2]);s.v[i].setZero();bulk.evaluate(i);}
   StdVec<Vecd>boundary={Vecd(-.05*dp,-.3*dp,t/2),Vecd(-.05*dp,.3*dp,t/2),Vecd(.05*dp,-.3*dp,t/2),Vecd(.05*dp,.3*dp,t/2)};
   StdVec<size_t>wrong_owner;for(const auto&q:boundary)wrong_owner.push_back(anchors({q},q[0]>0?-1:1)[0]);
   VisibleCrackSurface owned(port,boundary,wrong_owner,ref,dp,3,2,{},{},{},VisibleCrackSurface::AnchorOwnership::ReferencePoint);
   VisibleCrackSurface fixed_owner(port,boundary,wrong_owner,ref,dp,3,2);
   StdVec<Vecd>before,target;for(size_t q=0;q<boundary.size();++q){before.push_back(owned.position(q));target.push_back(before.back()+Vecd(0,0,.01*dp));}
   Real spring_area=dp*t/boundary.size();
   auto spring=[&](bool forces){Real U=0;if(forces)port.clear();for(size_t q=0;q<boundary.size();++q){Vecd d=owned.position(q)-target[q];U+=.5*K*spring_area*d.squaredNorm();if(forces)owned.scatter(q,-K*spring_area*d);}return U;};
   Real old_energy=bulk.energy()+spring(false),old_hg=bulk.stabilizationEnergy();
   owned.captureTransferValues();fixed_owner.captureTransferValues();
   ref.setTransferState(StdVec<Vecd>(s.x,s.x+count),StdVec<Matd>(bulk.F,bulk.F+count),StdVec<Vecd>(s.v,s.v+count));
   size_t actual_cuts=ref.cutPlane(Vecd::Zero(),Vecd::UnitX());bulk.rebuild();owned.rebuildChanged();
   bool fixed_rejected=must_reject([&](){fixed_owner.rebuildChanged();});
   Real transfer_error=0;size_t cross_entries=0;
   for(size_t q=0;q<boundary.size();++q){transfer_error=std::max(transfer_error,(owned.position(q)-before[q]).norm()/dp);if(boundary[q][0]*s.x0[owned.materialAnchor(q)][0]<=0)++cross_entries;for(size_t j:owned.rows[q].index)if(boundary[q][0]*s.x0[j][0]<=0||!ref.visibleReference(boundary[q],s.x0[j]))++cross_entries;}
   for(size_t i=0;i<count;++i)bulk.evaluate(i);
   Real energy_error=std::abs(bulk.energy()+spring(false)-old_energy)/std::max(std::abs(old_energy),Real(1e-30)),hg_loss=old_hg-bulk.stabilizationEnergy();
   VisibleCrackSurface restored(port,boundary,wrong_owner,ref,dp,3,2,{},owned.transferValues(),{},VisibleCrackSurface::AnchorOwnership::ReferencePoint);
   Real restart_error=0;for(size_t q=0;q<boundary.size();++q){if(restored.materialAnchor(q)!=owned.materialAnchor(q)||restored.rows[q].index!=owned.rows[q].index)throw std::runtime_error("Boundary reference ownership did not restore deterministically");for(size_t a=0;a<owned.rows[q].weight.size();++a)restart_error=std::max(restart_error,std::abs(restored.rows[q].weight[a]-owned.rows[q].weight[a]));}
   spring(true);Real power=0;StdVec<Vecd>directions(count);for(size_t i=0;i<count;++i){directions[i]=Vecd(.012*std::sin(5*s.x0[i][0]/dp),.008*std::cos(3*s.x0[i][1]/dp),.004);power+=s.mass[i]*port.external[i].dot(directions[i]);}
   Real eps=1e-9;for(size_t i=0;i<count;++i)s.x[i]+=eps*directions[i];Real ep=spring(false);for(size_t i=0;i<count;++i)s.x[i]-=2*eps*directions[i];Real em=spring(false);for(size_t i=0;i<count;++i)s.x[i]+=eps*directions[i];
   Real work_error=std::abs((ep-em)/(2*eps)+power)/std::max(std::abs(power),Real(1e-20));
   StdVec<Vecd>held;for(size_t q=0;q<boundary.size();++q)held.push_back(owned.position(q));
   for(size_t i=0;i<count;++i)if(s.x0[i][0]<0)s.x[i]+=Vecd(-20e-6,7e-6,0);
   Real opposite_motion=0;for(size_t q=0;q<boundary.size();++q)if(boundary[q][0]>0)opposite_motion=std::max(opposite_motion,(owned.position(q)-held[q]).norm()/dp);
   bool passed=fixed_rejected&&actual_cuts>0&&owned.reassignedPoints()==boundary.size()&&cross_entries==0&&transfer_error<1e-9&&energy_error<1e-8&&work_error<1e-5&&opposite_motion<1e-12&&restart_error<1e-12;
   std::ofstream report(out/"surface_ownership_acceptance.json");report<<std::setprecision(17)<<"{\"status\":\""<<(passed?"PASSED":"FAILED")<<"\",\"same_specimen_particles\":"<<count<<",\"reassigned_boundary_points\":"<<owned.reassignedPoints()<<",\"wrong_fixed_owner_rejected\":"<<(fixed_rejected?"true":"false")<<",\"cross_crack_interpolation_entries\":"<<cross_entries<<",\"position_transfer_error_over_spacing\":"<<transfer_error<<",\"physical_and_constraint_energy_relative_error\":"<<energy_error<<",\"separate_stabilizer_loss_J\":"<<hg_loss<<",\"force_work_relative_error\":"<<work_error<<",\"opposite_fragment_motion_leak_over_spacing\":"<<opposite_motion<<",\"restored_weight_error\":"<<restart_error<<",\"parameters_changed\":false}\n";report.close();
   for(size_t i=0;i<count;++i){s.x[i]=saved_x[i];s.v[i]=saved_v[i];}ref.restoreTopology(saved);bulk.rebuild();
   if(!passed)throw std::runtime_error("Same-specimen boundary ownership/energy acceptance failed");
  }
   auto automatic_growth=runAutomaticGrowthBenchmark(out,ref,bulk,port,intact_topology,s.x0,s.mass,s.x,s.v,acc,E,nu,ft,Gc,dp,t,K,nx,ny);
 if(automatic_growth.completed_only){std::ofstream report(out/"acceptance.json");report<<R"({"geometry_status":"PASSED","overall_status":"DYNAMIC_STAGE_COMPLETE_GC_PENDING","capsule_restart_allowed":false})"<<'\n';std::cout<<"COMPLETED_DYNAMIC_CHECKPOINT_SAVED"<<std::endl;return 86;}
  bool pass=geometry_pass&&automatic_growth.passed;
 std::ofstream report(out/"acceptance.json");report<<std::setprecision(17)<<"{\"geometry_status\":\""<<(geometry_pass?"PASSED":"FAILED")<<"\",\"overall_status\":\""<<(pass?"PASSED":"FAILED")<<"\",\"capsule_restart_allowed\":false,\"growth_insertion_status\":\""<<(automatic_growth.passed?"PASSED":"FAILED")<<"\",\"fracture_checkpoint_roundtrip_passed\":true,\"corrupt_or_incompatible_checkpoint_rejected\":true,\"topology_cache_guards_passed\":true,\"specimens\":1,\"physical_particles\":"<<count<<",\"disconnected_components\":"<<component_count<<",\"broken_bond_pairs\":"<<broken<<",\"live_damage_history_preserved\":"<<(damage_history_preserved?"true":"false")<<",\"live_damage_energy_over_E_volume\":"<<live_damage_energy<<",\"live_damage_force_over_E_area\":"<<live_damage_force<<",\"live_bonds_crossing_precrack\":"<<bad_bonds<<",\"maximum_visual_opening_m\":0.0006,\"sliding_m\":"<<2*dp<<",\"bulk_energy_over_E_volume\":"<<max_bulk<<",\"contact_force_error_over_E_area\":"<<max_contact_error<<",\"sliding_force_over_E_area\":"<<max_slip<<",\"gap_error_over_spacing\":"<<max_gap_error<<",\"force_energy_relative_error\":"<<max_fd<<",\"force_balance\":"<<max_balance<<",\"torque_balance\":"<<max_torque<<",\"outside_projections\":"<<outside<<",\"explicit_undamped_steps\":200,\"explicit_dt_s\":"<<dt<<",\"explicit_energy_relative_drift\":"<<drift<<",\"cohesive_Gc_J_m2\":"<<work<<",\"cohesive_Gc_relative_error\":"<<gc_error<<",\"cohesive_peak_relative_error\":"<<strength_error<<",\"growth_topology_insertion_tested\":true,\"material_parameters_changed\":false,\"particles_deleted\":0}\n";
 std::cout<<(pass?"PHYSICAL_PRECRACK_AND_AUTOMATIC_GROWTH_PASSED":"PHYSICAL_PRECRACK_BENCHMARK_FAILED")<<std::endl;return pass?0:2;
}catch(const std::exception&e){if(argc==2&&published){std::filesystem::path out(argv[1]);std::ifstream stage(out/"precrack_stage_acceptance.json");std::string status((std::istreambuf_iterator<char>(stage)),{});std::ofstream failed(out/"acceptance.json");failed<<"{\"overall_status\":\"FAILED\",\"geometry_status\":\""<<(status.find("PASSED")!=std::string::npos?"PASSED":"FAILED")<<"\",\"capsule_restart_allowed\":false,\"growth_insertion_status\":\"NOT_VERIFIED\",\"error\":"<<std::quoted(e.what())<<"}\n";}std::cerr<<"PHYSICAL_CRACK_BENCHMARK_ERROR "<<e.what()<<std::endl;return 1;}}
