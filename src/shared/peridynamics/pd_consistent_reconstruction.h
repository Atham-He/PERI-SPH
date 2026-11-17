#pragma once
#include "data_type.h"
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
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
    // Reorthogonalize small non-polynomial fields without subtracting two
    // large coordinate fields again. Normal-equation cancellation must not
    // turn an exact transfer constraint into apparent missing support.
    for(int pass=0;pass<2;++pass)residual-=B*(inverse*(B.transpose()*W.asDiagonal()*residual));
    Eigen::VectorXd scale=Eigen::VectorXd::Zero(Y.cols());
    for(int j=0;j<Y.cols();++j){Real data_scale=Y.col(j).norm()/std::sqrt(Real(Y.rows()));Real rms=std::sqrt((residual.col(j).array().square()*W.array()).sum());if(rms>1024*std::numeric_limits<Real>::epsilon()*std::max(data_scale,Real(1e-30)))scale[j]=1/rms;}
    M Z=residual*scale.asDiagonal(),weighted=W.cwiseSqrt().asDiagonal()*Z;
    // Direct singular values avoid squaring the conditioning of transferred
    // fields that can be almost polynomial or mutually dependent.
    Eigen::JacobiSVD<M> svd(weighted,Eigen::ComputeThinU|Eigen::ComputeThinV);svd.setThreshold(1e-12);
    M right_inverse=W.cwiseSqrt().asDiagonal()*svd.matrixU()*svd.singularValues().unaryExpr([&](Real v){return v>1e-12*svd.singularValues().maxCoeff()?1/v:0;}).asDiagonal()*svd.matrixV().transpose()*scale.asDiagonal();
    M next=G;
    for(int iteration=0;iteration<5;++iteration){M remaining=target-Y.transpose()*next;if(remaining.norm()<=tolerance)break;M correction=right_inverse*remaining;
        // Iterative refinement of the SAME constrained linear problem.
        correction-=W.asDiagonal()*B*(inverse*(B.transpose()*correction));
        if(!correction.allFinite())throw std::runtime_error("Nonfinite conservative reconstruction correction");next+=correction;
    }
    M change=next-G;Real polynomial_error=(B.transpose()*change).norm(),transfer_error=(Y.transpose()*next-target).norm();
    if(!next.allFinite()||polynomial_error>1e-8*std::max(Real(1),G.norm())||transfer_error>tolerance){
        if(const char*path=std::getenv("PERI_TRANSFER_AUDIT_FILE")){std::ifstream existing(path);if(!existing.good()){std::ofstream f(path);f<<std::setprecision(17)<<B.rows()<<' '<<B.cols()<<' '<<Y.cols()<<' '<<G.cols()<<' '<<tolerance<<'\n';auto matrix=[&](const auto&a){f<<a.rows()<<' '<<a.cols()<<'\n'<<a<<'\n';};matrix(B);matrix(W);matrix(Y);matrix(G);matrix(target);matrix(next);matrix(svd.singularValues());}}
        std::ostringstream message;message<<std::setprecision(17)<<"Conservative reconstruction unresolved: state_error="<<transfer_error<<", tolerance="<<tolerance<<", polynomial_error="<<polynomial_error<<", fields="<<Y.cols()<<", rank="<<svd.rank();throw std::runtime_error(message.str());
    }
    G=std::move(next);
}
}
