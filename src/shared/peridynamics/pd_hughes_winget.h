#pragma once

#include <Eigen/Dense>

namespace SPH::peridynamics
{
template <class Matrix>
EIGEN_DEVICE_FUNC Matrix hughesWingetStress(const Matrix &gradient, const Matrix &old,
                                          typename Matrix::Scalar lambda, typename Matrix::Scalar shear)
{
    const Matrix symmetric = (gradient + gradient.transpose()) * 0.5;
    const Matrix skew = (gradient - gradient.transpose()) * 0.5;
    const Matrix increment = lambda * symmetric.trace() * Matrix::Identity() + 2.0 * shear * symmetric;
    const Matrix rotation_system = Matrix::Identity() - skew * 0.5;
    const Matrix rotation = Matrix::Identity() + rotation_system.inverse() * skew;
    const Matrix rotated = rotation * old * rotation.transpose();
    return increment + rotated;
}
}
