#pragma once
#include "pd_extrinsic_crack_growth.h"
#include "pd_fracture_topology_checkpoint.h"
struct AutomaticGrowthAcceptance {bool passed=false;size_t new_facets=0;Real max_gap=0,energy_error=0,accounted_error=0,Gc_error=0;};
inline AutomaticGrowthAcceptance runAutomaticGrowthBenchmark(
    const std::filesystem::path&out,PDLayeredFractureReference&ref,PDBrokenBondEnergy&bulk,
    SolidForcePort&port,const PDLayeredFractureReference::TopologyState&intact,
    const StdVec<Vecd>&X,const Real*mass,Vecd*position,Vecd*velocity,Vecd*acc,
    Real E,Real nu,Real ft,Real Gc,Real dp,Real thickness,Real K,int nx,int ny)
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
    Growth growth(ref,bulk,port,ft,Gc,dp,K,bank);
    Field force(count),v(count,Vecd::Zero()),inverse_mass(count),old_force;
    std::vector<bool>fixed(count);for(size_t i=0;i<count;++i){size_t ix=(i%columns)/ny;fixed[i]=ix==0||ix==size_t(nx-1);inverse_mass[i]=fixed[i]?Vecd::Zero():Vecd::Constant(1/mass[i]);}
    auto bulkForce=[&](){for(size_t i=0;i<count;++i)bulk.evaluate(i);port.clear();growth.applyContact();for(size_t i=0;i<count;++i){acc[i].setZero();bulk.forces(i);force[i]=mass[i]*acc[i]+mass[i]*port.external[i];}};
    auto energy=[&](){bulkForce();growth.statistics(growth.positions());return bulk.energy()+growth.stored+growth.contact_energy;};
    auto refresh=[&](){growth.rebuildTraces();};
    energy();if(!growth.insert(0,true,refresh,energy))throw std::runtime_error("Automatic benchmark precrack: "+growth.last_rejection);
    Real initial=energy(),initialD=growth.dissipated;
    Real rho=1200,c=std::sqrt(E*(1-nu)/((1+nu)*(1-2*nu)*rho));
    Real ramp=8*L/c,end=ramp+2*L/c,maximum_displacement=.08*L;
    auto load=[&](Real time){Real a=std::clamp(time/ramp,Real(0),Real(1));return .5*maximum_displacement*(1-std::cos(Pi*a));};
    auto loadVelocity=[&](Real time){if(time>=ramp)return Real(0);return .5*maximum_displacement*Pi/ramp*std::sin(Pi*time/ramp);};
    Real elapsed=0,W=0,algorithmic=0,cohesive_force_work=0,max_error=0,max_accounted=0,peakW=0;size_t steps=0;bool fracture_restart_passed=false;
    std::ofstream csv(out/"automatic_growth.csv"),events(out/"automatic_growth_events.csv");csv<<std::setprecision(17)<<"step,time_s,grip_displacement_m,reaction_N,new_facets,failed_quadrature,max_gap_m,bulk_J,cohesive_J,contact_J,fracture_work_J,kinetic_J,boundary_work_J,cohesive_algorithmic_loss_J,cohesive_force_work_J,stabilizer_transfer_loss_J,insertion_defect_J,raw_budget_J,accounted_budget_J\n";events<<std::setprecision(17)<<"step,time_s,candidate,strength_ratio,new_area_m2,insertion_defect_J\n";
    std::filesystem::create_directories(out/"growth_frames");std::ofstream pvd(out/"automatic_growth.pvd");pvd<<"<VTKFile type=\"Collection\" version=\"0.1\"><Collection>\n";
    auto snapshot=[&](const std::string&name){StdVec<Vecd>q(position,position+count);StdVec<Real>damage(count,0);for(size_t i=0;i<count;++i){size_t dead=0;for(size_t n=0;n<ref.rows[i].size();++n)dead+=!ref.live(i,n);damage[i]=Real(dead)/ref.rows[i].size();}writeCloud(out/"growth_frames"/(name+".vtp"),q,damage);pvd<<"<DataSet timestep=\""<<elapsed<<"\" file=\"growth_frames/"<<name<<".vtp\"/>\n";};
    snapshot("initial");
    while(elapsed<end){
        if(steps%50==0){size_t k=growth.select();if(k<bank.size()){Real ratio=growth.strengthRatio(k),before=growth.insertion_defect,A=0;for(Real a:bank[k].area)A+=a;
            if(!growth.insert(k,false,refresh,energy)){auto a=growth.last_insertion;std::ofstream diagnostic(out/"insertion_rejected.csv");diagnostic<<std::setprecision(17)<<"candidate,area,old_energy,new_energy,old_bulk,new_bulk,old_cohesive,new_cohesive,old_contact,new_contact,new_dissipation,birth_gap,physical_defect,fraction,old_HG,new_HG,total_defect\n"<<a.candidate<<','<<a.area<<','<<a.old_energy<<','<<a.new_energy<<','<<a.old_bulk<<','<<a.new_bulk<<','<<a.old_cohesive<<','<<a.new_cohesive<<','<<a.old_contact<<','<<a.new_contact<<','<<a.new_dissipation<<','<<a.maximum_birth_gap<<','<<a.defect<<','<<a.fraction<<','<<a.old_stabilizer<<','<<a.new_stabilizer<<','<<a.total_defect<<'\n';throw std::runtime_error("Automatic benchmark insertion rejected: "+growth.last_rejection);}
            events<<steps<<','<<elapsed<<','<<k<<','<<ratio<<','<<A<<','<<growth.insertion_defect-before<<'\n';events.flush();snapshot("insert_"+std::to_string(k));
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
    result.passed=contact_cycle_passed&&fracture_restart_passed&&Gc_force_work_error<.01&&result.new_facets==3&&result.max_gap>2*growth.delta_f&&Gc_error<1e-6&&max_error<.02&&max_accounted<.005&&growth.max_insertion_fraction<.01;
    std::ofstream report(out/"automatic_growth_acceptance.json");report<<std::setprecision(17)<<"{\"status\":\""<<(result.passed?"PASSED":"FAILED")<<"\",\"same_specimen_particles\":"<<count<<",\"new_automatic_facets\":"<<result.new_facets<<",\"max_gap_m\":"<<growth.max_gap<<",\"contact_cycle_passed\":"<<(contact_cycle_passed?"true":"false")<<",\"contact_force_error\":"<<contact_error<<",\"contact_force_energy_error\":"<<contact_fd<<",\"contact_torque_error\":"<<contact_torque<<",\"postfracture_stressfree_bulk_over_EV\":"<<stressfree_bulk<<",\"neighbor_patch_contacts\":"<<neighbor_contacts<<",\"fracture_state_restart_passed\":"<<(fracture_restart_passed?"true":"false")<<",\"Gc_force_work_relative_error\":"<<Gc_force_work_error<<",\"Gc_work_relative_error\":"<<Gc_error<<",\"maximum_raw_energy_relative_error\":"<<max_error<<",\"maximum_accounted_energy_relative_error\":"<<max_accounted<<",\"maximum_insertion_energy_fraction\":"<<growth.max_insertion_fraction<<",\"steps\":"<<steps<<",\"damping\":0,\"parameter_changes\":false}\n";
    return result;
}
