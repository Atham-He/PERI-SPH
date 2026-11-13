/* Read-only observer. No force, material history, boundary, or dt feedback. */
#pragma once
#include "solid_snsni_plastic.h"
#include <filesystem>
#include <iomanip>
#include <sstream>
namespace SPH::solid_dynamics
{
class SnsniDiagnostics
{
 public:
    using Field=StdVec<Vecd>;
    struct Work {Real center=0,stabilization=0,nodal=0,absolute_scale=0;
        Real error()const{return std::abs(center+stabilization-nodal)/std::max(absolute_scale,Real(1e-30));}};
    static std::vector<Matd> centerMatrices(const SnsniReference&r,const SnsniPlasticSolid&b)
    {std::vector<Matd>A(r.cells.size());tbb::parallel_for(size_t(0),A.size(),[&](size_t l){A[l]=b.F[l].determinant()*r.cells[l].volume*b.pending[l].material.sigma*b.inverseF[l].transpose();});return A;}
    static Field centerForces(const SnsniReference&r,const SnsniPlasticSolid&b)
    {auto A=centerMatrices(r,b);Field f(r.X.size(),Vecd::Zero());tbb::parallel_for(size_t(0),f.size(),[&](size_t j){for(auto[l,k]:r.transpose[j])f[j]+=A[l]*r.cells[l].rows[k].gradient;});return f;}
    static Work instantVirtualWork(const SnsniReference&r,const SnsniPlasticSolid&b,const Field&v)
    {
        auto center=centerMatrices(r,b);std::vector<Work> local(r.cells.size());
        tbb::parallel_for(size_t(0),r.cells.size(),[&](size_t l){Matd g=Matd::Zero();std::array<Matd,3>h={Matd::Zero(),Matd::Zero(),Matd::Zero()};
            for(const auto&e:r.cells[l].rows){g+=v[e.j]*e.gradient.transpose();for(int a=0;a<3;++a)h[a]+=v[e.j][a]*e.hessian;}
            auto&w=local[l];w.center=(center[l].array()*g.array()).sum();w.stabilization=((b.force_gradient[l]-center[l]).array()*g.array()).sum();
            for(int a=0;a<3;++a)w.stabilization+=(b.force_hessian[l][a].array()*h[a].array()).sum();w.absolute_scale=std::abs(w.center)+std::abs(w.stabilization);});
        Work out;for(auto&w:local){out.center+=w.center;out.stabilization+=w.stabilization;out.absolute_scale+=w.absolute_scale;}
        for(size_t j=0;j<v.size();++j){Real w=b.internal[j].dot(v[j]);out.nodal+=w;out.absolute_scale+=std::abs(w);}return out;
    }
 private:
    struct Baseline {Real strain,upper,lower;};
    std::vector<Baseline>baseline;
    const SnsniReference&r;
    std::filesystem::path out;
    Real diameter,total_reference_volume=0,elapsed=0;
    Real old_kinetic=0,old_contact_energy=0,old_compression=0,old_mean=0;
    Field old_d,old_internal,old_center,old_total,center_cache;
    StdVec<Real>old_kappa;
    bool anchored=false,active=false;
    size_t matched=0;
    Real max_reaction_error=0,max_work_error=0,max_native_plastic_work_error=0;
    // Window cumulative values; only start after the accepted 8% state.
    Real Wcenter=0,Wstab=0,Wstab_positive=0,Wstab_negative=0,Wnodal=0,Wplate=0,Wexternal=0,Wdamping=0,delta_kinetic=0,delta_contact=0;
    std::array<Real,4>plastic_v{},plastic_q{},plastic_native{},branch_volume_time{},hardening_shear_work{};
    Real window_time=0;
    std::ofstream curve,budget,branches;
 public:
    SnsniDiagnostics(const std::filesystem::path&folder,const SnsniReference&ref,Real D,const std::string&reference_curve)
      :r(ref),out(folder/"diagnostics"),diameter(D)
    {
        std::filesystem::create_directories(out);
        std::ifstream f(reference_curve);if(!f)throw std::runtime_error("Missing original curve for diagnostic reproduction gate");
        std::string line;std::getline(f,line);std::vector<std::string>names;std::istringstream header(line);for(std::string x;std::getline(header,x,',');)names.push_back(x);
        auto col=[&](const std::string&name){auto i=std::find(names.begin(),names.end(),name);if(i==names.end())throw std::runtime_error("Baseline column missing");return size_t(i-names.begin());};
        size_t si=col("strain"),ui=col("upper_N"),li=col("lower_N");
        while(std::getline(f,line)){if(line.empty())continue;std::replace(line.begin(),line.end(),',',' ');std::istringstream input(line);std::vector<Real>v;for(Real x;input>>x;)v.push_back(x);baseline.push_back({v.at(si),v.at(ui),v.at(li)});}
        if(baseline.size()!=20)throw std::runtime_error("Expected all 20 original converged points");
        for(auto&c:r.cells)total_reference_volume+=c.volume;
        curve.open(out/"curve_reproduction.csv");curve<<std::setprecision(17)<<"stage,strain,baseline_upper_N,baseline_lower_N,diagnostic_upper_N,diagnostic_lower_N,max_absolute_difference_N,passed\n";
        budget.open(out/"work_budget.csv");budget<<std::setprecision(17)<<"step,stage,strain,elapsed_s,window_time_s,center_work_J,stabilization_work_J,stabilization_positive_J,stabilization_negative_J,internal_nodal_work_J,platen_work_J,external_particle_work_J,damping_loss_J,kinetic_change_J,contact_energy_change_J,budget_residual_J,contact_mapping_residual_J,integration_residual_J,max_split_relative_error,active_reference_area_top_m2,active_reference_area_bottom_m2,active_top_points,active_bottom_points,max_penetration_m,force_weighted_penetration_m,max_contact_pressure_Pa,contact_potential_J\n";
        branches.open(out/"branch_work.csv");branches<<std::setprecision(17)<<"step,stage,strain,branch,plastic_volume_work_J,plastic_deviatoric_work_J,native_plastic_work_J,volume_time_average,shear_Halpha_gamma_yield_component_J,max_native_work_difference_J\n";
        curve.flush();budget.flush();branches.flush();
        std::ofstream description(out/"definitions.txt");description<<"Read-only instrumentation; never feeds forces, histories, loading or dt.\nWork is trapezoidal generalized force times actual coefficient displacement increment.\nCenter force = pullback of J V0 sigma F^-T; stabilization = existing internal force minus center, independently checked using the gradient/Hessian weak form.\nNative branch work uses J V0 p delta(kappa) and J V0 q gamma, with signed volumetric work. Branch 0=elastic,1=shear,2=CAP,3=tensile apex.\nH alpha gamma is recorded only on the shear branch as a yield-work component, not stored hardening energy or a total hardening-force fraction.\nActive area is the reference quadrature area used in the penalty law. Snapshot CSV permits a separate deformed footprint-area calculation; do not confuse the two.\nAll work-budget residuals are reported without force correction.\n";
    }
    template<class Contact>void before(const SnsniCentralDifference&t,const SnsniPlasticSolid&b,const Contact&contact,const Field&total,size_t stage,Real compression)
    {
        active=anchored&&stage>=17;if(!active)return;
        old_d=t.displacement;old_internal=b.internal;old_total=total;
        old_center=center_cache;old_kinetic=t.kinetic();old_contact_energy=contact.energy;old_compression=compression;old_mean=.5*(contact.upper+contact.lower);
        old_kappa.resize(b.history.size());for(size_t l=0;l<b.history.size();++l)old_kappa[l]=b.history[l].kappa;
    }
    template<class Contact>void after(const SnsniCentralDifference&t,const SnsniPlasticSolid&b,const Contact&contact,const Field&total,Real dt,Real damping,size_t stage,Real compression)
    {
        elapsed+=dt;if(!active)return;
        center_cache=centerForces(r,b);Field increment(t.displacement.size());
        Real wc=0,ws=0,wi=0,we=0;
        for(size_t j=0;j<increment.size();++j){increment[j]=t.displacement[j]-old_d[j];
            wc+=.5*(old_center[j]+center_cache[j]).dot(increment[j]);
            ws+=.5*(old_internal[j]-old_center[j]+b.internal[j]-center_cache[j]).dot(increment[j]);
            wi+=.5*(old_internal[j]+b.internal[j]).dot(increment[j]);
            we+=.5*(old_total[j]+old_internal[j]+total[j]+b.internal[j]).dot(increment[j]);}
        Wcenter+=wc;Wstab+=ws;Wstab_positive+=std::max(ws,Real(0));Wstab_negative+=std::min(ws,Real(0));Wnodal+=wi;Wexternal+=we;
        Wplate+=.5*(old_mean+.5*(contact.upper+contact.lower))*(compression-old_compression);
        Real kinetic=t.kinetic(),factor=std::exp(-.5*damping*dt),factor2=factor*factor;
        Wdamping+=old_kinetic*(1-factor2)+kinetic*(1/factor2-1);delta_kinetic+=kinetic-old_kinetic;delta_contact+=contact.energy-old_contact_energy;
        for(size_t l=0;l<b.history.size();++l){const auto&m=b.pending[l].material;size_t branch=size_t(m.surface);if(branch>3)throw std::runtime_error("Unknown native DPC branch");
            Real JV=b.F[l].determinant()*r.cells[l].volume;
            Real wv=JV*m.p*(m.state.kappa-old_kappa[l]),wq=JV*m.q*m.plastic_shear;
            Real native=r.cells[l].volume*m.plastic_stress_work;
            max_native_plastic_work_error=std::max(max_native_plastic_work_error,std::abs(native-wv-wq));
            plastic_v[branch]+=wv;plastic_q[branch]+=wq;plastic_native[branch]+=native;branch_volume_time[branch]+=dt*r.cells[l].volume/total_reference_volume;
            if(branch==1)hardening_shear_work[branch]+=JV*b.model.shearHardening()*m.state.alpha*m.plastic_shear;
        }
        window_time+=dt;
        if(t.steps%200==0){Work check=instantVirtualWork(r,b,increment);max_work_error=std::max(max_work_error,check.error());
            if(check.error()>1e-10&&std::abs(check.center+check.stabilization-check.nodal)>1e-20)throw std::runtime_error("Diagnostic split failed independent weak-form check");
            writeBudget(t,b,contact,stage,compression);}
    }
 private:
    template<class Contact>void writeBudget(const SnsniCentralDifference&t,const SnsniPlasticSolid&,const Contact&c,size_t stage,Real compression)
    {
        Real atop=0,abottom=0,fwgap=0,sumf=0,maxgap=0,maxpressure=0;size_t ntop=0,nbottom=0;
        for(size_t q=0;q<c.force.size();++q){Real f=c.force[q][2],fn=std::abs(f);if(fn==0)continue;
            if(f<0){atop+=c.areas[q];++ntop;}else{abottom+=c.areas[q];++nbottom;}
            Real gap=fn/c.stiffness[q];fwgap+=fn*gap;sumf+=fn;maxgap=std::max(maxgap,gap);maxpressure=std::max(maxpressure,fn/c.areas[q]);}
        budget<<t.steps<<','<<stage<<','<<compression/diameter<<','<<elapsed<<','<<window_time<<','<<Wcenter<<','<<Wstab<<','<<Wstab_positive<<','<<Wstab_negative<<','<<Wnodal<<','<<Wplate<<','<<Wexternal<<','<<Wdamping<<','<<delta_kinetic<<','<<delta_contact<<','
          <<Wplate-Wnodal-delta_contact-delta_kinetic-Wdamping<<','<<Wplate-Wexternal-delta_contact<<','<<Wexternal-Wnodal-delta_kinetic-Wdamping<<','<<max_work_error<<','<<atop<<','<<abottom<<','<<ntop<<','<<nbottom<<','<<maxgap<<','<<(sumf>0?fwgap/sumf:0)<<','<<maxpressure<<','<<c.energy<<'\n';budget.flush();
        for(size_t branch=0;branch<4;++branch)branches<<t.steps<<','<<stage<<','<<compression/diameter<<','<<branch<<','<<plastic_v[branch]<<','<<plastic_q[branch]<<','<<plastic_native[branch]<<','<<branch_volume_time[branch]/window_time<<','<<hardening_shear_work[branch]<<','<<max_native_plastic_work_error<<'\n';branches.flush();
    }
    template<class Contact>void snapshot(size_t stage,const SnsniCentralDifference&t,const SnsniPlasticSolid&b,const Contact&c,Real stage_time,Real compression)
    {
        static_assert(sizeof(Real)==8,"Checkpoint requires double precision");
        std::string tag="stage_"+std::to_string(stage);std::ofstream data(out/(tag+"_state.bin"),std::ios::binary);
        auto raw=[&](const auto&x){data.write(reinterpret_cast<const char*>(&x),sizeof(x));};
        auto vector=[&](const Vecd&v){for(int a=0;a<3;++a)raw(v[a]);};
        auto matrix=[&](const Matd&m){for(int a=0;a<3;++a)for(int k=0;k<3;++k)raw(m(a,k));};
        uint64_t version=1,n=t.mass.size(),nl=b.history.size(),step=t.steps,load=stage;raw(version);raw(n);raw(nl);raw(step);raw(load);raw(elapsed);raw(stage_time);raw(compression);
        for(size_t j=0;j<n;++j){vector(t.displacement[j]);vector(t.velocity[j]);vector(t.acceleration[j]);raw(t.mass[j]);vector(b.internal[j]);}
        for(size_t l=0;l<nl;++l){const auto&m=b.pending[l].material;matrix(b.history[l].inverse_Cp);raw(b.history[l].kappa);raw(b.history[l].alpha);raw(b.accumulated_plastic_work[l]);matrix(b.F[l]);matrix(m.sigma);
            for(Real v:{m.p,m.q,m.relative_density,m.effective_d,m.plastic_shear,m.plastic_stress_work,m.elastic_energy,m.yield_residual})raw(v);
            int64_t surface=m.surface;raw(surface);matrix(b.force_gradient[l]);for(auto&h:b.force_hessian[l])matrix(h);matrix(b.correction[l]);matrix(b.pending[l].trial.sigma);for(auto&g:b.pending[l].trial.gradient)matrix(g);}
        if(!data)throw std::runtime_error("Checkpoint write failed");
        std::ofstream meta(out/(tag+"_state.json"));meta<<std::setprecision(17)<<"{\"format\":\"SNSNI_AUDIT_PRIMARY_STATE_V1\",\"endianness\":\"little\",\"stage\":"<<stage<<",\"step\":"<<step<<",\"elapsed_s\":"<<elapsed<<",\"stage_time_s\":"<<stage_time<<",\"compression_m\":"<<compression<<",\"next_stage\":"<<stage+1<<",\"next_stage_time_s\":0,\"coefficients\":"<<n<<",\"material_sites\":"<<nl<<",\"reader_restart_qualification\":false,\"per_coefficient_doubles\":13,\"per_material_doubles\":119,\"per_material_int64\":1}\n";
        // F/sigma/metrics and actual contact samples remain directly inspectable.
        std::ofstream material(out/(tag+"_material.csv"));material<<std::setprecision(17)<<"site,x0,y0,z0,x,y,z,J,p_Pa,q_Pa,kappa,alpha,relative_density,d_effective_Pa,pb_Pa,last_return_branch\n";
        for(size_t l=0;l<nl;++l){const auto&m=b.pending[l].material;auto a=b.model.parameters(m.relative_density,m.state.alpha);Vecd x=r.cells[l].center+r.value(l,t.displacement);
            material<<l<<','<<r.cells[l].center[0]<<','<<r.cells[l].center[1]<<','<<r.cells[l].center[2]<<','<<x[0]<<','<<x[1]<<','<<x[2]<<','<<b.F[l].determinant()<<','<<m.p<<','<<m.q<<','<<b.history[l].kappa<<','<<b.history[l].alpha<<','<<m.relative_density<<','<<m.effective_d<<','<<a.pb<<','<<m.surface<<'\n';}
        std::ofstream contacts(out/(tag+"_contact.csv"));contacts<<std::setprecision(17)<<"point,x0,y0,z0,x,y,z,fx_N,fy_N,fz_N,reference_area_m2,penetration_m\n";
        for(size_t q=0;q<c.force.size();++q)contacts<<q<<','<<c.surface.points[q][0]<<','<<c.surface.points[q][1]<<','<<c.surface.points[q][2]<<','<<c.position[q][0]<<','<<c.position[q][1]<<','<<c.position[q][2]<<','<<c.force[q][0]<<','<<c.force[q][1]<<','<<c.force[q][2]<<','<<c.areas[q]<<','<<std::abs(c.force[q][2])/c.stiffness[q]<<'\n';
        auto provenance=out.parent_path().parent_path()/"audit_provenance.json";
        if(std::filesystem::exists(provenance))std::filesystem::copy_file(provenance,out/(tag+"_provenance.json"),std::filesystem::copy_options::overwrite_existing);
    }
 public:
    template<class Contact>void equilibrium(size_t stage,const SnsniCentralDifference&t,const SnsniPlasticSolid&b,const Contact&c,Real stage_time,Real compression)
    {
        if(stage<1||stage>baseline.size())throw std::runtime_error("Diagnostic stage mismatch");auto ref=baseline[stage-1];
        Real error=std::max(std::abs(c.upper-ref.upper),std::abs(c.lower-ref.lower));bool pass=error<=1e-5&&std::abs(compression/diameter-ref.strain)<1e-13;
        curve<<stage<<','<<compression/diameter<<','<<ref.upper<<','<<ref.lower<<','<<c.upper<<','<<c.lower<<','<<error<<','<<(pass?1:0)<<'\n';curve.flush();max_reaction_error=std::max(max_reaction_error,error);matched+=pass;
        std::ofstream record(out/"reproduction_status.json");record<<std::setprecision(17)<<"{\"status\":\""<<(!pass?"FAILED":matched==20?"ALL_20_POINTS_MATCHED":"PREFIX_MATCHED")<<"\",\"matched_points\":"<<matched<<",\"planned_points\":20,\"max_absolute_reaction_difference_N\":"<<max_reaction_error<<",\"force_tolerance_N\":0.00001,\"diagnostic_interpretation_allowed\":"<<(matched==20?"true":"false")<<",\"max_work_split_error\":"<<max_work_error<<",\"max_native_plastic_work_error_J\":"<<max_native_plastic_work_error<<"}\n";record.close();
        if(!pass)throw std::runtime_error("Diagnostic curve differs from the original; no retry or parameter change");
        if(stage==16){anchored=true;center_cache=centerForces(r,b);}
        if(stage==16||stage==18||stage==20)snapshot(stage,t,b,c,stage_time,compression);
    }
};
}
