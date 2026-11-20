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
EIGEN_DEVICE_FUNC inline PointFailure evaluateHenckyTrial(const PointModel &model,const Tensor &F,const State &old,
    const TensorGradient &gradientF,const TensorGradient &gradientInverseCp,const Eigen::Vector3d &gradientKappa,
    TrialField &out, const Spectrum *prepared_spectrum = nullptr)
{
    const double J=F.determinant();
    if (!(J>0))
        return PointFailure::Deformation;
    const double density=model.initialDensity()*std::exp(old.kappa);
    Parameters a;
    const auto failure=model.parameters(density,old.alpha,a);
    if (failure!=PointFailure::None)
        return failure;
    const size_t k=model.densityInterval(density);
    const auto &left=model.node(k-1), &right=model.node(k);
    const double dE=(right.E-left.E)/(right.r-left.r);
    const double dnu=(right.nu-left.nu)/(right.r-left.r);
    const double dK=dE/(3*(1-2*a.nu))+2*a.E*dnu/(3*std::pow(1-2*a.nu,2));
    const double dG=dE/(2*(1+a.nu))-a.E*dnu/(2*std::pow(1+a.nu,2));
    Spectrum local_spectrum;
    if (!prepared_spectrum)
    {
        computeElasticSpectrum(F,old,local_spectrum);
        prepared_spectrum=&local_spectrum;
    }
    const auto &eig=*prepared_spectrum;
    if(eig.info()!=Eigen::Success||eig.eigenvalues().minCoeff()<=1e-14)
        return PointFailure::ElasticMetric;
    const auto&U=eig.eigenvectors();const Eigen::Vector3d values=eig.eigenvalues();
    Tensor strain=.5*U*values.array().log().matrix().asDiagonal()*U.transpose();
    Tensor dev=strain-strain.trace()/3*Tensor::Identity();
    Tensor tau=a.K*strain.trace()*Tensor::Identity()+2*a.G*dev;
    Tensor divided;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {double mean=.5*(values[i]+values[j]),z=(values[i]-values[j])/(values[i]+values[j]);
     divided(i,j)=std::abs(z)<1e-4?(1+z*z/3+std::pow(z,4)/5)/mean:
         (std::log(values[i])-std::log(values[j]))/(values[i]-values[j]);}
    out.sigma=tau/J;
    for(int axis=0;axis<3;++axis)
    {Tensor dbe=gradientF[axis]*old.inverse_Cp*F.transpose()+F*gradientInverseCp[axis]*F.transpose()+F*old.inverse_Cp*gradientF[axis].transpose();
     Tensor de=.5*U*(divided.array()*(U.transpose()*dbe*U).array()).matrix()*U.transpose();
     double dr=density*gradientKappa[axis];
     Tensor dtau=a.K*de.trace()*Tensor::Identity()+2*a.G*(de-de.trace()/3*Tensor::Identity())+
         dr*(dK*strain.trace()*Tensor::Identity()+2*dG*dev);
     out.gradient[axis]=dtau/J-out.sigma*(F.inverse()*gradientF[axis]).trace();}
    return PointFailure::None;
}
inline TrialField henckyTrialField(const Model &model,const Tensor &F,const State &old,
    const TensorGradient &gradientF,const TensorGradient &gradientInverseCp,const Eigen::Vector3d &gradientKappa)
{
    TrialField out;
    const auto failure=evaluateHenckyTrial(PointModel(model.pointInput()),F,old,gradientF,gradientInverseCp,
                                          gradientKappa,out);
    if (failure==PointFailure::Deformation)
        throw std::runtime_error("Invalid trial deformation");
    if (failure==PointFailure::ElasticMetric)
        throw std::runtime_error("Invalid trial metric");
    if (failure!=PointFailure::None)
        throw std::runtime_error(Model::failureReason(failure));
    return out;
}

EIGEN_DEVICE_FUNC inline double acousticStep(const Tensor &F, const Tensor &inverse_F, const Parameters &a,
                                            double rho0, double cfl_length)
{
    Spectrum spectrum;
    computeSymmetricSpectrum(inverse_F.transpose()*inverse_F,spectrum);
    const double inverse_stretch=std::sqrt(spectrum.eigenvalues().maxCoeff());
    const double sound=std::sqrt((a.K+4*a.G/3)*F.determinant()/rho0);
    return cfl_length/(inverse_stretch*sound);
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
