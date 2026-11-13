// DPC adapter for the trial-minus-correction construction of Wang et al.,
// S-NSNI Eqs.66,72,74-76. It retains the existing multiplicative Hencky DPC
// return map; it does not substitute a Hughes-Winget or isochoric J2 model.
// The constitutive correction must subsequently be reconstructed at integration
// sites with the RK basis. This file alone is not a plastic S-NSNI integrator.
#pragma once
#include "drucker_prager_cap_kernel.h"
#include <array>
namespace SPH::dpc
{
using TensorGradient=std::array<Tensor,3>;
struct TrialField
{
    Tensor sigma=Tensor::Zero();
    TensorGradient gradient{Tensor::Zero(),Tensor::Zero(),Tensor::Zero()};
};
struct SnsniMaterialStep
{
    Result material;
    TrialField trial;
    Tensor correction=Tensor::Zero();
};
inline TrialField henckyTrialField(const Model&model,const Tensor&F,const State&old,
    const TensorGradient&gradientF,const TensorGradient&gradientInverseCp,const Eigen::Vector3d&gradientKappa)
{
    const double J=F.determinant();
    if(!(J>0))throw std::runtime_error("Invalid trial deformation");
    const double density=model.r0()*std::exp(old.kappa);
    const auto a=model.parameters(density,old.alpha);
    const auto&nodes=model.nodes();size_t k=1;while(k+1<nodes.size()&&nodes[k].r<density)++k;
    const double dE=(nodes[k].E-nodes[k-1].E)/(nodes[k].r-nodes[k-1].r);
    const double dnu=(nodes[k].nu-nodes[k-1].nu)/(nodes[k].r-nodes[k-1].r);
    const double dK=dE/(3*(1-2*a.nu))+2*a.E*dnu/(3*std::pow(1-2*a.nu,2));
    const double dG=dE/(2*(1+a.nu))-a.E*dnu/(2*std::pow(1+a.nu,2));
    Tensor be=F*old.inverse_Cp*F.transpose();
    Eigen::SelfAdjointEigenSolver<Tensor> eig(.5*(be+be.transpose()));
    if(eig.info()!=Eigen::Success||eig.eigenvalues().minCoeff()<=1e-14)throw std::runtime_error("Invalid trial metric");
    const auto&U=eig.eigenvectors();const Eigen::Vector3d values=eig.eigenvalues();
    Tensor strain=.5*U*values.array().log().matrix().asDiagonal()*U.transpose();
    Tensor dev=strain-strain.trace()/3*Tensor::Identity();
    Tensor tau=a.K*strain.trace()*Tensor::Identity()+2*a.G*dev;
    Tensor divided;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {double mean=.5*(values[i]+values[j]),z=(values[i]-values[j])/(values[i]+values[j]);
     divided(i,j)=std::abs(z)<1e-4?(1+z*z/3+std::pow(z,4)/5)/mean:
         (std::log(values[i])-std::log(values[j]))/(values[i]-values[j]);}
    TrialField out;out.sigma=tau/J;
    for(int axis=0;axis<3;++axis)
    {Tensor dbe=gradientF[axis]*old.inverse_Cp*F.transpose()+F*gradientInverseCp[axis]*F.transpose()+F*old.inverse_Cp*gradientF[axis].transpose();
     Tensor de=.5*U*(divided.array()*(U.transpose()*dbe*U).array()).matrix()*U.transpose();
     double dr=density*gradientKappa[axis];
     Tensor dtau=a.K*de.trace()*Tensor::Identity()+2*a.G*(de-de.trace()/3*Tensor::Identity())+
         dr*(dK*strain.trace()*Tensor::Identity()+2*dG*dev);
     out.gradient[axis]=dtau/J-out.sigma*(F.inverse()*gradientF[axis]).trace();}
    return out;
}
inline SnsniMaterialStep snsniMaterialStep(const Model&model,const Tensor&F,const State&old,
    const TensorGradient&gradientF,const TensorGradient&gradientInverseCp,const Eigen::Vector3d&gradientKappa)
{
    SnsniMaterialStep step;
    step.trial=henckyTrialField(model,F,old,gradientF,gradientInverseCp,gradientKappa);
    step.material=model.update(F,old); // exactly one native return-map call
    step.correction=step.trial.sigma-step.material.sigma;
    return step; // caller controls the single history commit
}
}
