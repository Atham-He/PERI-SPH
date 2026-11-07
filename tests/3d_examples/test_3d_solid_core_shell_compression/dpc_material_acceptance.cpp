#include "drucker_prager_cap_kernel.h"
#include <Eigen/Geometry>
#include <iomanip>
#include <iostream>
using namespace SPH::dpc;
void require(bool value,const char *message){if(!value)throw std::runtime_error(message);}
int main(int argc,char **argv)
{
    try {
        require(argc==2,"provide density CSV"); const Model model=Model::fromCSV(argv[1],.79);
        const auto initial=model.parameters(.79); State zero;
        auto elastic=model.update(std::exp(-1e-5/3)*Tensor::Identity(),zero);
        const double expected=initial.K*1e-5/std::exp(-1e-5);
        require(std::abs(elastic.p/expected-1)<1e-8 && elastic.surface==0 && std::abs(elastic.state.kappa)<1e-12,"elastic bulk response");
        auto hydro=[&](int steps){State s;Result r;for(int k=1;k<=steps;++k){r=model.update(std::exp(-.06*k/(3*steps))*Tensor::Identity(),s);s=r.state;}return r;};
        const auto cap=hydro(120),fine=hydro(240);
        require(cap.surface==2 && cap.state.kappa>0 && cap.relative_density>.79,"hydrostatic cap must compact");
        require(std::abs(cap.p/model.parameters(cap.relative_density).pb-1)<2e-8,"hydrostatic cap pressure");
        require(std::abs(cap.p/fine.p-1)<1e-7 && std::abs(cap.state.kappa-fine.state.kappa)<1e-8,"hydrostatic increment consistency");
        require(std::abs(std::log(cap.state.inverse_Cp.determinant())-2*cap.state.kappa)<1e-9,"plastic volume/metric identity");
        // Unload analytically to zero elastic volume: total J equals plastic J.
        const auto unloaded=model.update(std::exp(-cap.state.kappa/3)*Tensor::Identity(),cap.state);
        require(std::abs(unloaded.p)<1. && std::abs(unloaded.state.kappa-cap.state.kappa)<1e-10,"unloading preserves residual densification");
        Tensor F=Tensor::Identity();F(0,0)=std::exp(.003);F(1,1)=std::exp(-.003);
        const auto shear=model.update(F,zero);
        require(shear.surface==1 && shear.plastic_shear>0 && shear.state.kappa<0 && shear.plastic_stress_work>=0,"nonassociated shear and admissible work");
        const Tensor Q=Eigen::AngleAxisd(.71,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
        const auto rotated=model.update(Q*F,zero);
        require((rotated.sigma-Q*shear.sigma*Q.transpose()).norm()/shear.sigma.norm()<1e-9,"finite rotation objectivity");
        const auto repeated=model.update(F,shear.state);
        require((repeated.sigma-shear.sigma).norm()/shear.sigma.norm()<1e-7 && std::abs(repeated.state.kappa-shear.state.kappa)<1e-9,"stationary state must not flow repeatedly");
        bool domain_rejected=false;try{model.update(std::exp(-.5/3)*Tensor::Identity(),zero);}catch(const std::runtime_error&){domain_rejected=true;}
        require(domain_rejected,"no extrapolation past density table");
        std::cout<<std::setprecision(17)<<"{\"status\":\"PASSED_CONSTITUTIVE_ACCEPTANCE\",\"hydrostatic_p_Pa\":"<<cap.p
                 <<",\"hydrostatic_kappa\":"<<cap.state.kappa<<",\"residual_relative_density\":"<<cap.relative_density
                 <<",\"shear_kappa\":"<<shear.state.kappa<<",\"shear_plastic_stress_work_Pa\":"<<shear.plastic_stress_work
                 <<",\"scope\":\"material integration only; no spatial or experimental validation\"}\n";
        return 0;
    }catch(const std::exception&e){std::cerr<<"DPC ACCEPTANCE FAILED: "<<e.what()<<'\n';return 1;}
}
