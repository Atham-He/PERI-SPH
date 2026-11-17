#pragma once
#include "data_type.h"
#include <Eigen/Eigenvalues>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>
namespace SPH::peridynamics {
/** Minimum weighted change of a discrete reconstruction, subject to its old
 * polynomial moments AND a transferred continuous-state functional.
 * No q, velocity, mass, material constant or force is altered. The resulting
 * operator is frozen between topology events and is used with its exact
 * transpose; affine/quadratic patch consistency and frame objectivity remain.
 * This is a constrained GMLS state transfer, not an elastic-energy offset.
 *
 * B: nodal polynomial values; Y: three transferred coordinate-field values;
 * G: old new-support functional (one value column or three gradient columns).
 * Find dG minimizing ||W^-1/2 dG|| with B^T dG=0 and
 * Y^T(G+dG)=target. Polynomial residualization avoids dependent affine modes.
 */
inline void transferReconstruction(const Eigen::MatrixXd&B,const Eigen::VectorXd&W,
    const Eigen::MatrixXd&Y,Eigen::MatrixXd&G,const Eigen::MatrixXd&target,Real tolerance)
{
    using M=Eigen::MatrixXd;
    if(B.rows()!=W.size()||B.rows()!=Y.rows()||B.rows()!=G.rows()||target.rows()!=Y.cols()||target.cols()!=G.cols()||!B.allFinite()||!W.allFinite()||!Y.allFinite()||!G.allFinite()||!target.allFinite()||(W.array()<=0).any())throw std::invalid_argument("Invalid conservative reconstruction dimensions/data");
    M error=target-Y.transpose()*G;if(error.norm()<=tolerance)return;
    M moment=B.transpose()*W.asDiagonal()*B;
    Eigen::SelfAdjointEigenSolver<M> eig(moment);if(eig.info()!=Eigen::Success)throw std::runtime_error("Constrained GMLS polynomial factorization");
    auto values=eig.eigenvalues();M inverse=M::Zero(moment.rows(),moment.cols());Real largest=values.maxCoeff();
    for(int i=0;i<values.size();++i)if(values[i]>1e-12*largest)inverse.noalias()+=eig.eigenvectors().col(i)*eig.eigenvectors().col(i).transpose()/values[i];
    M residual=Y-B*(inverse*(B.transpose()*W.asDiagonal()*Y));
    Eigen::VectorXd scale=Eigen::VectorXd::Zero(Y.cols());
    for(int j=0;j<Y.cols();++j){Real data_scale=Y.col(j).norm()/std::sqrt(Real(Y.rows()));Real rms=std::sqrt((residual.col(j).array().square()*W.array()).sum());if(rms>1024*std::numeric_limits<Real>::epsilon()*std::max(data_scale,Real(1e-30)))scale[j]=1/rms;}
    M Z=residual*scale.asDiagonal();M gram=Z.transpose()*W.asDiagonal()*Z;
    Eigen::SelfAdjointEigenSolver<M> e(gram);if(e.info()!=Eigen::Success)throw std::runtime_error("Constrained GMLS transfer factorization");
    M pinv=M::Zero(Y.cols(),Y.cols());Real max=e.eigenvalues().maxCoeff();for(int j=0;j<Y.cols();++j)if(e.eigenvalues()[j]>1e-11*max)pinv.noalias()+=e.eigenvectors().col(j)*e.eigenvectors().col(j).transpose()/e.eigenvalues()[j];
    M correction=W.asDiagonal()*Z*pinv*scale.asDiagonal()*error;
    if(!correction.allFinite()||(B.transpose()*correction).norm()>1e-8*std::max(Real(1),G.norm()))throw std::runtime_error("Conservative reconstruction violates polynomial consistency");
    M next=G+correction;
    if((Y.transpose()*next-target).norm()>tolerance)throw std::runtime_error("Live support cannot transfer the continuous state consistently; no forces or positions were corrected");
    G=std::move(next);
}
}
