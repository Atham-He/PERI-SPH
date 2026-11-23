#pragma once

#include "pd_finite_volume_biot_batch.h"
#include "pd_finite_volume_biot_terms.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <limits>

namespace SPH::peridynamics::cuda_detail
{
__device__ inline bool finiteBiotTensor(const Eigen::Matrix3d &tensor)
{
    for (int i=0;i<9;++i)
        if (!isfinite(tensor.data()[i]))
            return false;
    return true;
}

__device__ inline FiniteVolumeBiotPointResult evaluateFiniteVolumeBiotPoint(const double (&input)[3][3],double bulk,double shear)
{
    FiniteVolumeBiotPointResult result{};
    Eigen::Matrix3d deformation;
    for (int a=0;a<3;++a)
        for (int b=0;b<3;++b)
        {
            const double value=input[a][b];
            if (!isfinite(value)) return result;
            deformation(a,b)=value;
        }
    const double J=deformation.determinant();
    if (!(J>0) || !isfinite(J)) return result;

    const Eigen::Matrix3d metric=deformation.transpose()*deformation;
    if (!finiteBiotTensor(metric)) return result;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> spectrum;
    spectrum.computeDirect(metric);
    const Eigen::Vector3d lambda=spectrum.eigenvalues();
    const double relative_error=128*std::numeric_limits<double>::epsilon();
    if (spectrum.info()!=Eigen::Success||!(lambda[0]>relative_error*lambda[2])) return result;
    const auto &axes=spectrum.eigenvectors();
    const Eigen::Matrix3d eigen_residual=metric*axes-axes*lambda.asDiagonal();
    const Eigen::Matrix3d orthogonality=axes.transpose()*axes-Eigen::Matrix3d::Identity();
    if (!(eigen_residual.cwiseAbs().maxCoeff()<=relative_error*metric.cwiseAbs().maxCoeff()) ||
        !(orthogonality.cwiseAbs().maxCoeff()<=relative_error))
        return result;

    const Eigen::Vector3d stretch=lambda.cwiseSqrt();
    const Eigen::Vector3d dev=stretch-Eigen::Vector3d::Constant(stretch.mean());
    const Eigen::Vector3d strain_over_stretch=dev.cwiseQuotient(stretch);
    const auto volume=finiteVolumeBiotVolumetric(bulk,J);
    // If F = U diag(stretch) V^T, this is U diag(dev) V^T.
    Eigen::Matrix3d stretched_axes=deformation*axes;
    for (int a=0;a<3;++a)
        for (int b=0;b<3;++b)
            stretched_axes(a,b)*=strain_over_stretch[b];
    const Eigen::Matrix3d distortional=stretched_axes*axes.transpose();
    const Eigen::Matrix3d inverse=deformation.inverse();
    const Eigen::Matrix3d piola=2*shear*distortional+volume.piola_scale*inverse.transpose();
    result.jacobian=J;
    result.elastic_energy=shear*dev.squaredNorm()+volume.energy;
    result.tangent_bound=finiteVolumeBiotTangent(bulk,shear,stretch,J,0.);
    if (!finiteBiotTensor(piola)||!isfinite(result.elastic_energy)||!isfinite(result.tangent_bound)) return result;
    for (int a=0;a<3;++a)
        for (int b=0;b<3;++b)
            result.piola.derivative[a][b]=piola(a,b);
    result.accepted=true;
    return result;
}
}
