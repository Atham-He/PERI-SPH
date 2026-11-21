#pragma once
#include "pd_extrinsic_crack_growth.h"
#include "pd_fracture_topology_checkpoint.h"
#include "solid_snsni_plastic.h"
#include "recovery_checkpoint.h"
struct AutomaticGrowthAcceptance {bool passed=false,completed_only=false;size_t new_facets=0;Real max_gap=0,energy_error=0,accounted_error=0,Gc_error=0;};
inline AutomaticGrowthAcceptance runAutomaticGrowthBenchmark(
    const std::filesystem::path&out,PDLayeredFractureReference&ref,PDBrokenBondEnergy&bulk,
    SolidForcePort&port,const PDLayeredFractureReference::TopologyState&intact,
    const StdVec<Vecd>&X,const Real*mass,Vecd*position,Vecd*velocity,Vecd*acc,
    Real E,Real nu,Real ft,Real Gc,Real dp,Real thickness,Real K,int nx,int ny,
    CrackProjectionBatch *projection, std::unique_ptr<CrackStrengthBatch> strength_backend)
{
    using Growth=ExtrinsicCrackGrowth;using Field=std::vector<Vecd>;
    const size_t count=X.size(),columns=count/3;const Real L=nx*dp,Y=ny*dp;
    ref.restoreTopology(intact);bulk.rebuild();for(size_t i=0;i<count;++i){position[i]=X[i];velocity[i].setZero();}
    std::vector<Growth::Candidate>bank;
    ThicknessQuadrature rule(3);
    for(size_t segment=0;segment<4;++segment){Growth::Candidate c;c.id=segment;Real lo=-Y/2+segment*Y/4,hi=lo+Y/4;c.normal=Vecd::UnitX();
        c.vertex={Vecd(0,lo,-thickness/2),Vecd(0,hi,-thickness/2),Vecd(0,hi,thickness/2),Vecd(0,lo,thickness/2)};
        c.barrier={{c.vertex[0],c.vertex[1],c.vertex[2]},{c.vertex[0],c.vertex[2],c.vertex[3]}};
        for(Real a:{-1/std::sqrt(Real(3)),1/std::sqrt(Real(3))})for(size_t k=0;k<3;++k){c.point.push_back(Vecd(0,.5*(hi+lo)+.5*(hi-lo)*a,.5*thickness*rule.xi[k]));c.area.push_back(.5*(hi-lo)*.5*thickness*rule.weight[k]);}
        for(size_t layer=0;layer<3;++layer){size_t minus=count,plus=count;Real dm=1,dp2=1;Vecd q(0,.5*(lo+hi),.5*thickness*rule.xi[layer]);for(size_t j=layer*columns;j<(layer+1)*columns;++j){Real distance=(X[j]-q).squaredNorm();if(X[j][0]<0&&distance<dm){dm=distance;minus=j;}if(X[j][0]>0&&distance<dp2){dp2=distance;plus=j;}}c.stress_pairs.push_back({minus,plus});}bank.push_back(c);
    }
    Growth growth(ref,bulk,port,ft,Gc,dp,K,bank,Growth::SurfaceLaw::Cohesive,projection,std::move(strength_backend));
    Field force(count),v(count,Vecd::Zero()),inverse_mass(count),old_force;
    std::vector<bool>fixed(count);for(size_t i=0;i<count;++i){size_t ix=(i%columns)/ny;fixed[i]=ix==0||ix==size_t(nx-1);inverse_mass[i]=fixed[i]?Vecd::Zero():Vecd::Constant(1/mass[i]);}
    auto bulkForce=[&](){bulk.evaluateAll();port.clear();growth.applyContact();bulk.recomputeForces();tbb::parallel_for(size_t(0),count,[&](size_t i){force[i]=mass[i]*acc[i]+mass[i]*port.external[i];});};
    auto energy=[&](){bulkForce();growth.statistics(growth.positions());return bulk.energy()+growth.stored+growth.contact_energy;};
    auto refresh=[&](){growth.rebuildTraces();};
    energy();if(!growth.insert(0,true,refresh,energy))throw std::runtime_error("Automatic benchmark precrack: "+growth.last_rejection);
    Real initial=energy(),initialD=growth.dissipated;
    Real rho=1200,c=std::sqrt(E*(1-nu)/((1+nu)*(1-2*nu)*rho));
    Real ramp=8*L/c,end=ramp+2*L/c,maximum_displacement=.08*L;
    auto load=[&](Real time){Real a=std::clamp(time/ramp,Real(0),Real(1));return .5*maximum_displacement*(1-std::cos(Pi*a));};
    auto loadVelocity=[&](Real time){if(time>=ramp)return Real(0);return .5*maximum_displacement*Pi/ramp*std::sin(Pi*time/ramp);};
    Real elapsed=0,W=0,algorithmic=0,cohesive_force_work=0,max_error=0,max_accounted=0,peakW=0;size_t steps=0;bool fracture_restart_passed=false;
    std::ofstream csv(out/"automatic_growth.csv"),events(out/"automatic_growth_events.csv");csv<<std::setprecision(17)<<"step,time_s,grip_displacement_m,reaction_N,new_facets,failed_quadrature,max_gap_m,bulk_J,cohesive_J,contact_J,fracture_work_J,kinetic_J,boundary_work_J,cohesive_algorithmic_loss_J,cohesive_force_work_J,stabilizer_transfer_loss_J,insertion_defect_J,raw_budget_J,accounted_budget_J\n";events<<std::setprecision(17)<<"step,time_s,candidate,strength_ratio,new_area_m2,insertion_defect_J,bulk_dt_s,stabilizer_dt_s,contact_dt_s,max_gradient_norm_per_m,max_trace_l1\n";
    bool replay_active=std::getenv("PERI_GROWTH_REPLAY_STATE")!=nullptr;
    if(replay_active){std::ifstream frozen(std::getenv("PERI_GROWTH_REPLAY_STATE"));std::string line;std::getline(frozen,line);size_t read=0;while(std::getline(frozen,line)){std::replace(line.begin(),line.end(),',',' ');std::istringstream values(line);size_t id;values>>id;if(id!=read||id>=count)throw std::runtime_error("Invalid frozen diagnostic state index");for(int d=0;d<3;++d)values>>position[id][d];for(int d=0;d<3;++d)values>>v[id][d];if(!values)throw std::runtime_error("Invalid frozen diagnostic field");velocity[id]=v[id];++read;}if(read!=count)throw std::runtime_error("Incomplete diagnostic state");elapsed=std::stod(std::getenv("PERI_GROWTH_REPLAY_TIME"));steps=std::stoull(std::getenv("PERI_GROWTH_REPLAY_STEP"));Real kinetic=0;for(size_t i=0;i<count;++i)if(!fixed[i])kinetic+=.5*mass[i]*v[i].squaredNorm();W=energy()-initial+kinetic;peakW=std::abs(W);std::cout<<"DIAGNOSTIC_REPLAY_ONLY_NO_FORMAL_ACCEPTANCE"<<std::endl;}
    std::filesystem::create_directories(out/"growth_frames");std::ofstream pvd(out/"automatic_growth.pvd");pvd<<"<VTKFile type=\"Collection\" version=\"0.1\"><Collection>\n";
    auto snapshot=[&](const std::string&name){StdVec<Vecd>q(position,position+count);StdVec<Real>damage(count,0);for(size_t i=0;i<count;++i){size_t dead=0;for(size_t n=0;n<ref.rows[i].size();++n)dead+=!ref.live(i,n);damage[i]=Real(dead)/ref.rows[i].size();}writeCloud(out/"growth_frames"/(name+".vtp"),q,damage);pvd<<"<DataSet timestep=\""<<elapsed<<"\" file=\"growth_frames/"<<name<<".vtp\"/>\n";};
    // A complete dynamic-end state, not a trace-only codec: exact q/v,
    // operator history, cohesive history and work counters can resume the
    // prescribed-opening/contact checks without repeating the dynamic solve.
    std::string physics_identity=std::getenv("PERI_GROWTH_PHYSICS_FINGERPRINT")?std::getenv("PERI_GROWTH_PHYSICS_FINGERPRINT"):"self_roundtrip_only";
    auto completedTransfer=[&](recovery::Archive&a){
        std::string magic="NOSB_SLAB_COMPLETED_DYNAMIC_V1",identity=physics_identity;uint64_t schema=1;
        std::vector<Real>parameters{E,nu,ft,Gc,dp,thickness,K,Real(nx),Real(ny),ramp,end,maximum_displacement},expected_parameters=parameters;
        a.fields(magic,schema,identity,parameters);if(magic!="NOSB_SLAB_COMPLETED_DYNAMIC_V1"||schema!=1||identity!=physics_identity||parameters!=expected_parameters)throw std::runtime_error("Completed dynamic checkpoint identity mismatch");
        a.fields(initial,initialD,elapsed,W,algorithmic,cohesive_force_work,max_error,max_accounted,peakW,steps,fracture_restart_passed,growth.insertion_defect,growth.max_insertion_fraction,growth.stabilizer_transfer_loss,growth.insertion_events,growth.rejected_events,growth.initial_dissipation);
        Field q(position,position+count),velocity_state=v,traction=growth.last_traction;std::vector<size_t>order;StdVec<StdVec<Vecd>>minus_values,plus_values;
        if(!a.reading)for(const auto&patch:growth.patches){order.push_back(patch.candidate);minus_values.push_back(patch.minus->transferValues());plus_values.push_back(patch.plus->transferValues());}
        a.fields(q,velocity_state,traction,order,minus_values,plus_values);
        if(q.size()!=count||velocity_state.size()!=count)throw std::runtime_error("Completed state node count mismatch");for(size_t i=0;i<count;++i)if(!q[i].allFinite()||!velocity_state[i].allFinite())throw std::runtime_error("Nonfinite completed state field");
        using Archive=FractureTopologyCheckpoint;Archive::Histories histories;Archive::IMEXBlock block;block.name="complete_dynamic_state";block.ft=ft;block.Gc=Gc;
        for(size_t k:order){if(k>=bank.size())throw std::runtime_error("Completed state candidate outside bank");for(size_t j=0;j<bank[k].point.size();++j)block.history_ids.push_back(64*bank[k].id+j);}
        if(a.reading)block.history.resize(block.history_ids.size());else block.history=growth.history;histories.imex.push_back(block);std::vector<uint64_t>ids(count);for(size_t i=0;i<count;++i)ids[i]=i;
        Archive::Bytes bytes;if(!a.reading)bytes=Archive::encode(ref,physics_identity,ids,histories);a.value(bytes);
        if(a.reading){for(size_t i=0;i<count;++i){position[i]=q[i];velocity[i]=velocity_state[i];}v=velocity_state;Archive::load(ref,bytes,physics_identity,ids,histories);growth.restorePatches(order,histories.imex[0].history,minus_values,plus_values);if(traction.size()!=growth.history.size())throw std::runtime_error("Completed state traction count mismatch");growth.last_traction=traction;bulk.rebuild();}
        bulkForce();auto cohesive=growth.cohesiveForces(growth.positions());Field expected_force=force;for(size_t i=0;i<count;++i)expected_force[i]+=cohesive[i];Real expected_energy=energy();a.fields(expected_force,expected_energy);if(expected_force.size()!=count||!std::isfinite(expected_energy))throw std::runtime_error("Invalid completed state force/energy");for(const auto&f:expected_force)if(!f.allFinite())throw std::runtime_error("Nonfinite completed state force");
        if(a.reading){Real error=0,scale=0;for(size_t i=0;i<count;++i){error=std::max(error,(force[i]+cohesive[i]-expected_force[i]).norm());scale=std::max(scale,expected_force[i].norm());}Real U=energy();if(error>1e-9*std::max(scale,Real(1e-20))||std::abs(U-expected_energy)>1e-10*std::max(std::abs(expected_energy),Real(1e-30)))throw std::runtime_error("Completed state force/energy restore mismatch");}
        if(std::abs(elapsed-end)>1e-18)throw std::runtime_error("Checkpoint is not the completed dynamic stage");
    };
    auto loadCompleted=[&](const std::filesystem::path&file){recovery::Archive a(file,true);completedTransfer(a);a.finish();};
    bool completed_state_restart_passed=false;
    if(const char*file=std::getenv("PERI_GROWTH_RESUME_COMPLETED_DYNAMIC")){if(!std::getenv("PERI_GROWTH_PHYSICS_FINGERPRINT"))throw std::runtime_error("External completed-state resume requires a pinned physics fingerprint");loadCompleted(file);completed_state_restart_passed=true;std::cout<<"RESUMED_COMPLETED_DYNAMIC_STATE steps="<<steps<<std::endl;}
    snapshot(completed_state_restart_passed?"resumed_completed":"initial");size_t diagnostic_stop=replay_active?steps+(std::getenv("PERI_GROWTH_REPLAY_STEPS")?std::stoull(std::getenv("PERI_GROWTH_REPLAY_STEPS")):500):~size_t(0);
    while(elapsed<end){
        if(steps>=diagnostic_stop)throw std::runtime_error("Frozen-state diagnostic replay reached its requested limit; no formal acceptance issued");
        if(steps%50==0){for(size_t event=0;event<bank.size();++event){size_t k=growth.select();if(k==bank.size())break;Real ratio=growth.strengthRatio(k),before=growth.insertion_defect,A=0;{std::ofstream frozen(out/("preinsert_"+std::to_string(k)+".csv"));frozen<<std::setprecision(17)<<"node,x,y,z,vx,vy,vz\n";for(size_t j=0;j<count;++j)frozen<<j<<','<<position[j][0]<<','<<position[j][1]<<','<<position[j][2]<<','<<velocity[j][0]<<','<<velocity[j][1]<<','<<velocity[j][2]<<'\n';}for(Real a:bank[k].area)A+=a;
            if(!growth.insert(k,false,refresh,energy)){auto a=growth.last_insertion;std::ofstream diagnostic(out/"insertion_rejected.csv");diagnostic<<std::setprecision(17)<<"candidate,area,old_energy,new_energy,old_bulk,new_bulk,old_cohesive,new_cohesive,old_contact,new_contact,new_dissipation,birth_gap,physical_defect,fraction,old_HG,new_HG,total_defect\n"<<a.candidate<<','<<a.area<<','<<a.old_energy<<','<<a.new_energy<<','<<a.old_bulk<<','<<a.new_bulk<<','<<a.old_cohesive<<','<<a.new_cohesive<<','<<a.old_contact<<','<<a.new_contact<<','<<a.new_dissipation<<','<<a.maximum_birth_gap<<','<<a.defect<<','<<a.fraction<<','<<a.old_stabilizer<<','<<a.new_stabilizer<<','<<a.total_defect<<'\n';throw std::runtime_error("Automatic benchmark insertion rejected: "+growth.last_rejection);}
            Real l1=0;for(const auto&patch:growth.patches)for(const auto*surface:{patch.minus.get(),patch.plus.get()})for(size_t q=0;q<4+surface->natural.size();++q){Real value=0;for(Real w:surface->rows[q].weight)value+=std::abs(w);l1=std::max(l1,value);}events<<steps<<','<<elapsed<<','<<k<<','<<ratio<<','<<A<<','<<growth.insertion_defect-before<<','<<bulk.finiteStrainTimeStep()<<','<<bulk.stabilizer_dt<<','<<growth.contact_dt<<','<<ref.maximum_gradient_norm<<','<<l1<<'\n';events.flush();snapshot("insert_"+std::to_string(k));
            if(!fracture_restart_passed){
                using Archive=FractureTopologyCheckpoint;Archive::Histories histories;Archive::IMEXBlock block;block.name="automatic_ft_Gc_surface";block.ft=ft;block.Gc=Gc;block.history=growth.history;
                std::vector<size_t>order;StdVec<StdVec<Vecd>>minus_values,plus_values,minus_velocity,plus_velocity;
                for(const auto&p:growth.patches){order.push_back(p.candidate);minus_values.push_back(p.minus->transferValues());plus_values.push_back(p.plus->transferValues());minus_velocity.push_back(p.minus->transferVelocities());plus_velocity.push_back(p.plus->transferVelocities());for(size_t j=0;j<bank[p.candidate].point.size();++j)block.history_ids.push_back(64*bank[p.candidate].id+j);}
                histories.imex.push_back(block);std::vector<uint64_t>ids(count);for(size_t i=0;i<count;++i)ids[i]=i;
                auto bytes=Archive::encode(ref,"same_automatic_growth_state_transfer_v3",ids,histories);auto traction=growth.last_traction;bulkForce();Field old=force;auto old_cohesive=growth.cohesiveForces(growth.positions());for(size_t i=0;i<count;++i)old[i]+=old_cohesive[i];
                Archive::load(ref,bytes,"same_automatic_growth_state_transfer_v3",ids,histories);growth.restorePatches(order,histories.imex[0].history,minus_values,plus_values,minus_velocity,plus_velocity);growth.last_traction=traction;bulk.rebuild();bulkForce();auto restored=growth.cohesiveForces(growth.positions());
                Real error=0,scale=0;for(size_t i=0;i<count;++i){error=std::max(error,(force[i]+restored[i]-old[i]).norm());scale=std::max(scale,old[i].norm());}
                if(error>1e-9*std::max(scale,Real(1e-20))||bytes!=Archive::encode(ref,"same_automatic_growth_state_transfer_v3",ids,histories))throw std::runtime_error("Automatic fracture-state restore changed its force/operator history");
                std::ofstream saved(out/"automatic_fracture_state.bin",std::ios::binary);saved.write(reinterpret_cast<const char*>(bytes.data()),bytes.size());fracture_restart_passed=true;
            }
        }}
        bulkForce();Real bound=std::min({bulk.finiteStrainTimeStep(),bulk.stabilizer_dt,growth.contact_dt});Real dt=std::min(.25*bound,end-elapsed);
        if(!std::isfinite(dt)||dt<1e-16||steps>1000000)throw std::runtime_error("Automatic benchmark time bound exhausted");
        Field qold(position,position+count),vhalf=v,qfree=qold;old_force=force;Real psiold=growth.stored;
        for(size_t i=0;i<count;++i){if(!fixed[i]){vhalf[i]+=.5*dt*force[i]/mass[i];qfree[i]+=dt*vhalf[i];}else{qfree[i]=X[i]+Vecd((X[i][0]>0?.5:-.5)*load(elapsed+dt),0,0);}}
        auto solution=growth.drift(qfree,qold,inverse_mass,dt);
        Real velocity_defect=0,surface_work=0;for(size_t i=0;i<count;++i){position[i]=solution.positions[i];surface_work-=solution.force[i].dot(solution.positions[i]-qold[i]);if(!fixed[i]){Vecd after=vhalf[i]+dt*solution.force[i]/mass[i];velocity_defect+=.5*mass[i]*(after-vhalf[i]).squaredNorm();v[i]=after;}}
        cohesive_force_work+=surface_work;algorithmic+=velocity_defect+surface_work-(solution.stored-psiold)-solution.dissipation_increment;
        bulkForce();Real kinetic=0;for(size_t i=0;i<count;++i){if(!fixed[i]){v[i]+=.5*dt*force[i]/mass[i];kinetic+=.5*mass[i]*v[i].squaredNorm();}else{v[i]=Vecd((X[i][0]>0?.5:-.5)*loadVelocity(elapsed+dt),0,0);W-=.5*(old_force[i]+force[i]+solution.force[i]).dot(solution.positions[i]-qold[i]);}velocity[i]=v[i];}
        elapsed+=dt;++steps;growth.statistics(growth.positions());Real U=bulk.energy()+growth.stored+growth.contact_energy,D=growth.dissipated-initialD,raw=W-(U-initial)-D-kinetic,accounted=raw-algorithmic-growth.stabilizer_transfer_loss+growth.insertion_defect;peakW=std::max(peakW,std::abs(W));
        if(peakW>ft*ft/E*L*Y*thickness*.01){max_error=std::max(max_error,std::abs(raw)/peakW);max_accounted=std::max(max_accounted,std::abs(accounted)/peakW);}
        if(steps%200==0||elapsed>=end){Real reaction=0;for(size_t i=0;i<count;++i)if(fixed[i]&&X[i][0]>0)reaction-=force[i][0]+solution.force[i][0];csv<<steps<<','<<elapsed<<','<<load(elapsed)<<','<<reaction<<','<<growth.patches.size()-1<<','<<growth.failed_quadrature<<','<<growth.max_gap<<','<<bulk.energy()<<','<<growth.stored<<','<<growth.contact_energy<<','<<D<<','<<kinetic<<','<<W<<','<<algorithmic<<','<<cohesive_force_work<<','<<growth.stabilizer_transfer_loss<<','<<growth.insertion_defect<<','<<raw<<','<<accounted<<'\n';csv.flush();}
        if(steps%5000==0){std::cout<<"GROWTH step="<<steps<<" time="<<elapsed<<" facets="<<growth.patches.size()<<" gap="<<growth.max_gap<<" energy_error="<<max_error<<std::endl;snapshot("step_"+std::to_string(steps));}
    }
    snapshot("final");pvd<<"</Collection></VTKFile>\n";
    auto completed_file=out/"automatic_completed_state.bin";{auto temp=std::filesystem::path(completed_file.string()+".tmp");recovery::Archive a(temp,false);completedTransfer(a);a.finish();std::filesystem::rename(temp,completed_file);}loadCompleted(completed_file);completed_state_restart_passed=true;
    if(std::getenv("PERI_GROWTH_STOP_AFTER_COMPLETED_DYNAMIC")){
        if(replay_active||growth.patches.size()!=4||!fracture_restart_passed||max_error>=.02||max_accounted>=.005||growth.max_insertion_fraction>=.01)throw std::runtime_error("Completed dynamic stage failed its unchanged acceptance gates");
        csv.flush();events.flush();pvd.flush();AutomaticGrowthAcceptance partial;partial.completed_only=true;partial.new_facets=growth.patches.size()-1;partial.energy_error=max_error;partial.accounted_error=max_accounted;return partial;
    }
    auto quadratureAudit=[&](const std::filesystem::path&file){std::ofstream f(file);f<<std::setprecision(17)<<"candidate,quadrature,area_m2,kappa_m,normal_gap_m,effective_gap_m,nx,ny,nz\n";CohesiveSurfaceEvolution model(ft,Gc,growth.quadrature());auto q=growth.positions();size_t index=0;for(const auto&patch:growth.patches)for(size_t k=0;k<bank[patch.candidate].point.size();++k,++index){auto g=model.geometry(index,q);f<<patch.candidate<<','<<k<<','<<bank[patch.candidate].area[k]<<','<<growth.history[index].kappa<<','<<g.normal_gap<<','<<g.jump.norm()<<','<<g.normal.transpose().format(Eigen::IOFormat(17,0,",",",","","","",""))<<'\n';}};
    quadratureAudit(out/"dynamic_end_quadrature.csv");Real dynamic_end_gap=growth.max_gap,dynamic_end_fracture_work=growth.dissipated-initialD;
    // A free tensile strip can retain a subcritical ligament after unloading.
    // The SAME history therefore completes a prescribed rigid-fragment opening
    // for the independent spatial force-work/Gc check. This is a load-parameter
    // virtual-work test, not extra physical time or an extension of the dynamic
    // energy ledger. No history is set to failure: kappa advances only from the
    // actual separation. Physical masses and the saved dynamic velocities stay.
    Real controlled_translation=0,controlled_work=0,controlled_cohesive_work=0,controlled_energy_error=0;
    if(growth.failed_quadrature!=growth.history.size()){
        for(size_t i=0;i<count;++i)for(size_t k=0;k<ref.rows[i].size();++k)if(ref.live(i,k)&&X[i][0]*X[ref.rows[i][k].j][0]<0)throw std::runtime_error("Rigid-fragment opening requires actual completed material cuts");
        CohesiveSurfaceEvolution model(ft,Gc,growth.quadrature());auto base=growth.positions();Vecd direction=Vecd::Zero();for(size_t k=0;k<growth.history.size();++k)if(growth.history[k].kappa<growth.delta_f)direction+=model.quadrature[k].area*model.geometry(k,base).normal;
        if(direction.norm()==0)throw std::runtime_error("Unresolved prescribed opening direction");direction.normalize();
        for(size_t k=0;k<growth.history.size();++k)if(growth.history[k].kappa<growth.delta_f&&direction.dot(model.geometry(k,base).normal)<=0)throw std::runtime_error("Pending crack normals do not admit one physical rigid opening direction");
        controlled_translation=2*growth.delta_f;Real Ubegin=energy(),Dbegin=growth.dissipated,area=0;for(const auto&patch:growth.patches)if(patch.candidate)for(Real value:bank[patch.candidate].area)area+=value;
        const size_t increments=4000;std::ofstream table(out/"controlled_opening.csv");table<<std::setprecision(17)<<"increment,rigid_relative_opening_m,failed_quadrature,cohesive_work_J,total_constraint_work_J,recoverable_J,fracture_work_increment_J,potential_work_residual_J\n";
        std::filesystem::create_directories(out/"controlled_opening_frames");std::ofstream path(out/"controlled_opening.pvd");path<<"<VTKFile type=\"Collection\" version=\"0.1\"><Collection>\n";
        for(size_t step=1;step<=increments;++step){bulkForce();Field old_contact_bulk=force,old_cohesive=growth.cohesiveForces(growth.positions()),before=growth.positions();Real shift=controlled_translation*Real(step)/increments;
            for(size_t i=0;i<count;++i)position[i]=base[i]+(X[i][0]>0?.5:-.5)*shift*direction;
            auto q=growth.positions();for(size_t k=0;k<growth.history.size();++k)growth.history[k].kappa=std::max(growth.history[k].kappa,model.geometry(k,q).jump.norm());
            Real U=energy(),D=growth.dissipated;Field current_cohesive=growth.cohesiveForces(q);for(size_t i=0;i<count;++i){Vecd delta=q[i]-before[i];controlled_cohesive_work-=.5*(old_cohesive[i]+current_cohesive[i]).dot(delta);controlled_work-=.5*(old_contact_bulk[i]+force[i]+old_cohesive[i]+current_cohesive[i]).dot(delta);}
            Real defect=controlled_work-(U-Ubegin)-(D-Dbegin);controlled_energy_error=std::max(controlled_energy_error,std::abs(defect)/(Gc*area));
            if(step%100==0||step==increments)table<<step<<','<<shift<<','<<growth.failed_quadrature<<','<<controlled_cohesive_work<<','<<controlled_work<<','<<U<<','<<D-Dbegin<<','<<defect<<'\n';
            if(step%400==0||step==1){std::string tag=std::to_string(step);StdVec<Real>marker(count,0);for(size_t i=0;i<count;++i){size_t dead=0;for(size_t k=0;k<ref.rows[i].size();++k)dead+=!ref.live(i,k);marker[i]=Real(dead)/ref.rows[i].size();}writeCloud(out/"controlled_opening_frames"/("material_"+tag+".vtp"),q,marker);
                std::ofstream frame(out/"controlled_opening_frames"/("frame_"+tag+".vtm"));frame<<R"(<VTKFile type="vtkMultiBlockDataSet" version="1.0"><vtkMultiBlockDataSet>)";frame<<"<DataSet index=\"0\" name=\"material\" file=\"material_"<<tag<<".vtp\"/>";size_t part=1;
                for(const auto&patch:growth.patches)for(auto side:{std::pair{patch.minus.get(),"minus"},std::pair{patch.plus.get(),"plus"}}){StdVec<Vecd>corners;for(size_t k=0;k<4;++k)corners.push_back(side.first->position(k));std::string name="face_"+std::to_string(patch.candidate)+"_"+side.second+"_"+tag+".vtp";writeFace(out/"controlled_opening_frames"/name,corners);frame<<"<DataSet index=\""<<part++<<"\" name=\"crack_"<<patch.candidate<<'_'<<side.second<<"\" file=\""<<name<<"\"/>";}
                frame<<"</vtkMultiBlockDataSet></VTKFile>\n";path<<"<DataSet timestep=\""<<shift<<"\" file=\"controlled_opening_frames/frame_"<<tag<<".vtm\"/>\n";}
        }
        path<<"</Collection></VTKFile>\n";cohesive_force_work+=controlled_cohesive_work;
    }
    quadratureAudit(out/"controlled_opening_end_quadrature.csv");
    // The SAME fully fractured specimen now checks reopening, closure and
    // sliding across neighbouring patches. No damage history is reset/forced.
    bool contact_cycle_passed=false;Real contact_fd=0,contact_error=0,contact_torque=0,stressfree_bulk=0;size_t neighbor_contacts=0;
    Field saved_position(position,position+count),saved_velocity(velocity,velocity+count);
    if(growth.failed_quadrature==growth.history.size()&&growth.patches.size()==4){
        struct Mode{Real gap,slide,angle;};std::vector<Mode>modes={{20e-6,0,0},{-2e-9,0,0},{-2e-9,2*dp,0},{-2e-9,2*dp,.37},{20e-6,2*dp,.37}};
        Real volume=0;for(size_t i=0;i<count;++i)volume+=mass[i]/rho;
        for(const auto&m:modes){Matd Q=Eigen::AngleAxis<Real>(m.angle,Vecd(1,2,3).normalized()).toRotationMatrix();Field direction(count);
            for(size_t i=0;i<count;++i){Vecd shift((X[i][0]>0?.5:-.5)*m.gap,X[i][0]>0?m.slide:0,0);position[i]=Q*(X[i]+shift);direction[i]=Vecd(.01*std::sin(X[i][0]/dp),.008*std::cos(X[i][1]/dp),.004*std::sin(X[i][2]/thickness));velocity[i]=direction[i];}
            bulkForce();stressfree_bulk=std::max(stressfree_bulk,std::abs(bulk.energy())/(E*volume));auto&network=*growth.contact_network;Real expected=K*(Y-std::abs(m.slide))*thickness*std::max(-m.gap,Real(0));contact_error=std::max(contact_error,std::abs(network.normal_force_sum-expected)/std::max(K*Y*thickness*2e-9,Real(1e-30)));contact_torque=std::max(contact_torque,network.torque_residual);neighbor_contacts+=network.neighbor_contacts;
            if(m.gap<0){Real power=0;for(size_t i=0;i<count;++i)power+=force[i].dot(direction[i]);Real eps=1e-10;Field base(position,position+count);for(size_t i=0;i<count;++i)position[i]=base[i]+eps*direction[i];bulkForce();Real ep=bulk.energy()+growth.contact_energy;for(size_t i=0;i<count;++i)position[i]=base[i]-eps*direction[i];bulkForce();Real em=bulk.energy()+growth.contact_energy;contact_fd=std::max(contact_fd,std::abs((ep-em)/(2*eps)+power)/std::max(std::abs(power),Real(1e-12)));}
        }
        contact_cycle_passed=contact_error<1e-7&&contact_fd<1e-5&&contact_torque<1e-9&&stressfree_bulk<1e-12&&neighbor_contacts>0;
    }
    for(size_t i=0;i<count;++i){position[i]=saved_position[i];velocity[i]=saved_velocity[i];}bulkForce();growth.statistics(growth.positions());
    Real new_area=0;for(size_t k=1;k<bank.size();++k)if(growth.active[k])for(Real a:bank[k].area)new_area+=a;
    Real newD=growth.dissipated-initialD,Gc_force_work_error=new_area>0?std::abs(cohesive_force_work/(Gc*new_area)-1):1,Gc_error=new_area>0?std::abs(newD/(Gc*new_area)-1):1;
    AutomaticGrowthAcceptance result;result.new_facets=growth.patches.size()-1;result.max_gap=growth.max_gap;result.energy_error=max_error;result.accounted_error=max_accounted;result.Gc_error=Gc_error;
    result.passed=!replay_active&&completed_state_restart_passed&&controlled_energy_error<.005&&contact_cycle_passed&&fracture_restart_passed&&Gc_force_work_error<.01&&result.new_facets==3&&result.max_gap>2*growth.delta_f&&Gc_error<1e-6&&max_error<.02&&max_accounted<.005&&growth.max_insertion_fraction<.01;
    std::ofstream protocol(out/"opening_protocol.json");protocol<<std::setprecision(17)<<R"({"dynamic_history_unchanged":true,"dynamic_end_gap_m":)"<<dynamic_end_gap<<R"(,"dynamic_end_fracture_work_J":)"<<dynamic_end_fracture_work<<R"(,"completed_dynamic_state_restart_passed":)"<<(completed_state_restart_passed?"true":"false")<<R"(,"controlled_rigid_translation_m":)"<<controlled_translation<<R"(,"controlled_cohesive_work_J":)"<<controlled_cohesive_work<<R"(,"controlled_potential_work_relative_error":)"<<controlled_energy_error<<R"(,"controlled_opening_is_virtual_work_parameter_not_dynamic_time":true,"damage_set_by_separation_only":true,"physical_parameters_changed":false})"<<'\n';
    std::ofstream report(out/"automatic_growth_acceptance.json");report<<std::setprecision(17)<<"{\"status\":\""<<(result.passed?"PASSED":"FAILED")<<"\",\"same_specimen_particles\":"<<count<<",\"new_automatic_facets\":"<<result.new_facets<<",\"max_gap_m\":"<<growth.max_gap<<",\"contact_cycle_passed\":"<<(contact_cycle_passed?"true":"false")<<",\"contact_force_error\":"<<contact_error<<",\"contact_force_energy_error\":"<<contact_fd<<",\"contact_torque_error\":"<<contact_torque<<",\"postfracture_stressfree_bulk_over_EV\":"<<stressfree_bulk<<",\"neighbor_patch_contacts\":"<<neighbor_contacts<<",\"fracture_state_restart_passed\":"<<(fracture_restart_passed?"true":"false")<<",\"Gc_force_work_relative_error\":"<<Gc_force_work_error<<",\"Gc_work_relative_error\":"<<Gc_error<<",\"maximum_raw_energy_relative_error\":"<<max_error<<",\"maximum_accounted_energy_relative_error\":"<<max_accounted<<",\"maximum_insertion_energy_fraction\":"<<growth.max_insertion_fraction<<",\"steps\":"<<steps<<",\"damping\":0,\"parameter_changes\":false}\n";
    return result;
}
