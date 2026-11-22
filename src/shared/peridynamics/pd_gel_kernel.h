#pragma once

#include "pd_gel_batch.h"
#include <Eigen/Dense>
#include <cmath>

namespace SPH::peridynamics
{
template<class Matrix>
EIGEN_DEVICE_FUNC bool finiteGelTensor(const Matrix &tensor)
{
    using std::isfinite;
    for (int a = 0; a < Matrix::RowsAtCompileTime; ++a)
        for (int b = 0; b < Matrix::ColsAtCompileTime; ++b)
            if (!isfinite(tensor(a,b)))
                return false;
    return true;
}

template<class Matrix>
EIGEN_DEVICE_FUNC bool validGelDeformation(const Matrix &F, typename Matrix::Scalar J)
{
    using std::isfinite;
    return finiteGelTensor(F) && isfinite(J) && J >= 0.05;
}

template<class Matrix>
EIGEN_DEVICE_FUNC Matrix deviatoricGelGreen(const Matrix &F)
{
    const Matrix E = 0.5 * (F.transpose() * F - Matrix::Identity());
    return E - (E.trace() / Matrix::RowsAtCompileTime) * Matrix::Identity();
}

template<class Matrix>
EIGEN_DEVICE_FUNC bool validGelHistoryTensor(const Matrix &tensor)
{
    using std::abs;
    return finiteGelTensor(tensor) &&
           !((tensor.transpose() - tensor).norm() > 1e-12 * (1 + tensor.norm())) &&
           !(abs(tensor.trace()) > 1e-12 * (1 + tensor.norm()));
}

template<class History>
EIGEN_DEVICE_FUNC GelKernelStatus gelHistoryStatus(const History &history)
{
    using std::isfinite;
    if (!validGelHistoryTensor(history.green) || !validGelHistoryTensor(history.viscous))
        return GelKernelStatus::InvalidHistory;
    if (!isfinite(history.dissipation) || history.dissipation < 0)
        return GelKernelStatus::InvalidDissipation;
    return GelKernelStatus::OK;
}

template<class Matrix, class Parameters, class History>
EIGEN_DEVICE_FUNC GelKernelStatus advanceGelKernel(const Matrix &F, typename Matrix::Scalar dt,
                                                    const Parameters &p, const History &old, History &next)
{
    using Scalar = typename Matrix::Scalar;
    using std::exp;
    using std::expm1;
    using std::isfinite;
    const auto status = gelHistoryStatus(old);
    if (status != GelKernelStatus::OK)
        return status;
    if (!validGelDeformation(F, F.determinant()))
        return GelKernelStatus::InvalidDeformation;
    const Matrix green = deviatoricGelGreen(F);
    if (!isfinite(dt) || dt < 0)
        return GelKernelStatus::InvalidIncrement;
    if (dt == 0)
    {
        next = old;
        return GelKernelStatus::OK;
    }
    const Scalar z = dt / p.relaxation_time;
    if (!isfinite(z) || z <= 0)
        return GelKernelStatus::InvalidRelaxation;
    const Scalar decay = exp(-z), one_minus = -expm1(-z);
    const Matrix q0 = old.green - old.viscous, change = green - old.green;
    const Matrix q1 = decay * q0 + (one_minus / z) * change;
    // Integral of (1-exp(-s))^2 / z^2. The series avoids cancellation.
    const Scalar w22 = z < 1e-3 ? z * (1.0 / 3 + z * (-1.0 / 4 + z * (7.0 / 60 +
        z * (-1.0 / 24 + z * (31.0 / 2520 - z / 320))))) :
        (1 - (one_minus + 0.5 * one_minus * one_minus) / z) / z;
    const Scalar increment = 2 * p.branch_shear * (
        (-0.5 * expm1(-2 * z)) * q0.squaredNorm() +
        (one_minus * one_minus / z) * q0.cwiseProduct(change).sum() + w22 * change.squaredNorm());
    if (!isfinite(increment) || increment < 0)
        return GelKernelStatus::InvalidDissipationIntegral;
    next = {green, green - q1, old.dissipation + increment};
    return gelHistoryStatus(next);
}

template<class Matrix, class Parameters, class History, class Response>
EIGEN_DEVICE_FUNC GelKernelStatus gelResponseKernel(const Matrix &F, const Parameters &p,
                                                     const History &history, Response &result)
{
    using Scalar = typename Matrix::Scalar;
    using std::isfinite;
    using std::log;
    using std::pow;
    const auto status = gelHistoryStatus(history);
    if (status != GelKernelStatus::OK)
        return status;
    const Scalar J = F.determinant();
    if (!validGelDeformation(F, J))
        return GelKernelStatus::InvalidDeformation;
    const Scalar isochoric = pow(J, -2.0 / 3.0);
    const Scalar trace_C = F.squaredNorm();
    const Matrix inverse_transpose = F.inverse().transpose();
    const Matrix q = deviatoricGelGreen(F) - history.viscous;
    result.free_energy = 0.5 * p.equilibrium_shear * (isochoric * trace_C - 3) +
        0.25 * p.bulk * (J * J - 1 - 2 * log(J)) + p.branch_shear * q.squaredNorm();
    result.pk1 = p.equilibrium_shear * isochoric * (F - trace_C / 3 * inverse_transpose) +
        0.5 * p.bulk * (J * J - 1) * inverse_transpose + 2 * p.branch_shear * F * q;
    result.cauchy = result.pk1 * F.transpose() / J;
    if (!finiteGelTensor(result.pk1) || !finiteGelTensor(result.cauchy) || !isfinite(result.free_energy))
        return GelKernelStatus::NonfiniteResponse;
    return GelKernelStatus::OK;
}
}
