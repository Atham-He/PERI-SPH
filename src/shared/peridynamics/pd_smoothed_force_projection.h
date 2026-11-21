#pragma once

#include <Eigen/Core>
#include <Eigen/LU>
#include <array>

namespace SPH::peridynamics
{
using SmoothedForceTensor = Eigen::Matrix3d;
using SmoothedForceGradient = std::array<SmoothedForceTensor,3>;

// Contract the spatial stress gradient against the pushed-forward cell
// moment. The correction gradient is supplied in reference coordinates.
EIGEN_DEVICE_FUNC inline void projectSmoothedForce(
    const SmoothedForceTensor &F, const SmoothedForceTensor &inverse_F,
    double volume, const SmoothedForceTensor &reference_moment,
    const SmoothedForceTensor &stress, const SmoothedForceGradient &trial_gradient,
    const SmoothedForceGradient &reference_correction_gradient,
    const SmoothedForceGradient &displacement_hessian,
    SmoothedForceTensor &nodal_gradient, SmoothedForceGradient &nodal_hessian)
{
    SmoothedForceGradient sigma_gradient;
    for (int k=0;k<3;++k)
    {
        sigma_gradient[k]=trial_gradient[k];
        for (int j=0;j<3;++j)
            sigma_gradient[k]-=inverse_F(j,k)*reference_correction_gradient[j];
    }

    const double J=F.determinant();
    const SmoothedForceTensor moment=J*F*reference_moment*F.transpose();
    SmoothedForceTensor geometric=SmoothedForceTensor::Zero();
    for (int a=0;a<3;++a)
    {
        auto &coefficient=nodal_hessian[a];
        coefficient.setZero();
        for (int k=0;k<3;++k)
            coefficient+=(inverse_F*sigma_gradient[k].row(a).transpose())*
                         (inverse_F*moment.row(k).transpose()).transpose();
        for (int b=0;b<3;++b)
            geometric(a,b)=(coefficient.array()*displacement_hessian[b].array()).sum();
    }
    // The Hessian transformation retains the geometric term from the
    // displacement Hessian as well as the material stress correction.
    nodal_gradient=(J*volume*stress-geometric)*inverse_F.transpose();
}
}
