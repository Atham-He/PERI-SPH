#pragma once
/** Directional effective crack energy, following the QR mode relaxation of
 * Hakimzadeh et al., JMPS 167 (2022) 104994, Definition 2.1.
 * Adaptation to the EXISTING finite-volume Biot energy; E and nu are unchanged.
 * Crack-normal opening and frictionless sliding relax; tangential material
 * deformation and normal closure retain energy. No residual shear stiffness,
 * contact penalty, particle deletion or fracture-parameter rescaling is used.
 */
#include <array>
#include <cmath>
#include <limits>

namespace SPH::solid_dynamics::directional_detail
{
/** Exact first/second derivatives of the local energy in nine F components. */
struct Jet
{
    Real v;std::array<Real,9>g{};std::array<Real,81>h{};
    Jet(Real value=0):v(value){}
    static Jet variable(Real value,size_t i){Jet a(value);a.g[i]=1;return a;}
};
inline Jet operator+(const Jet&a,const Jet&b){Jet r(a.v+b.v);for(size_t i=0;i<9;++i)r.g[i]=a.g[i]+b.g[i];for(size_t i=0;i<81;++i)r.h[i]=a.h[i]+b.h[i];return r;}
inline Jet operator-(const Jet&a,const Jet&b){Jet r(a.v-b.v);for(size_t i=0;i<9;++i)r.g[i]=a.g[i]-b.g[i];for(size_t i=0;i<81;++i)r.h[i]=a.h[i]-b.h[i];return r;}
inline Jet operator*(const Jet&a,const Jet&b){Jet r(a.v*b.v);for(size_t i=0;i<9;++i)r.g[i]=a.g[i]*b.v+a.v*b.g[i];for(size_t i=0;i<9;++i)for(size_t j=0;j<9;++j)r.h[9*i+j]=a.h[9*i+j]*b.v+a.g[i]*b.g[j]+a.g[j]*b.g[i]+a.v*b.h[9*i+j];return r;}
inline Jet compose(const Jet&a,Real value,Real first,Real second){Jet r(value);for(size_t i=0;i<9;++i)r.g[i]=first*a.g[i];for(size_t i=0;i<9;++i)for(size_t j=0;j<9;++j)r.h[9*i+j]=first*a.h[9*i+j]+second*a.g[i]*a.g[j];return r;}
inline Jet inverse(const Jet&a){return compose(a,1/a.v,-1/(a.v*a.v),2/(a.v*a.v*a.v));}
inline Jet operator/(const Jet&a,const Jet&b){return a*inverse(b);}
inline Real value(Real a){return a;}inline Real value(const Jet&a){return a.v;}
inline Real root(Real a){return std::sqrt(a);}inline Jet root(const Jet&a){Real r=std::sqrt(a.v);return compose(a,r,.5/r,-.25/(a.v*r));}
inline Real logarithm1p(Real a){return std::log1p(a);}inline Jet logarithm1p(const Jet&a){return compose(a,std::log1p(a.v),1/(1+a.v),-1/((1+a.v)*(1+a.v)));}
template<class T>T bulk(const T&J,Real K)
{
    T z=J-T(1);
    // Stable evaluation of J^2-1-2 log(J), not a constitutive approximation.
    if(std::abs(value(z))<1e-3){T p=T(Real(1)/5);p=T(-Real(2)/9)+z*p;p=T(Real(1)/4)+z*p;p=T(-Real(2)/7)+z*p;p=T(Real(1)/3)+z*p;p=T(-Real(2)/5)+z*p;p=T(.5)+z*p;p=T(-Real(2)/3)+z*p;p=T(2)+z*p;return T(.25*K)*z*z*p;}
    return T(.25*K)*(z*(J+T(1))-T(2)*logarithm1p(z));
}
template<class T>struct Evaluation{T energy,normal_height,elastic_height;bool closed;};
template<class T>Evaluation<T> evaluate(const std::array<T,9>&F,const Vecd&N,Real K,Real mu)
{
    Vecd t=N.unitOrthogonal(),s=N.cross(t);std::array<T,3>u{},v{},w{},c{};
    for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j){u[i]=u[i]+F[3*i+j]*T(t[j]);v[i]=v[i]+F[3*i+j]*T(s[j]);w[i]=w[i]+F[3*i+j]*T(N[j]);}
    c[0]=u[1]*v[2]-u[2]*v[1];c[1]=u[2]*v[0]-u[0]*v[2];c[2]=u[0]*v[1]-u[1]*v[0];
    T uu(0),vv(0),uv(0),area2(0),det(0);for(size_t i=0;i<3;++i){uu=uu+u[i]*u[i];vv=vv+v[i]*v[i];uv=uv+u[i]*v[i];area2=area2+c[i]*c[i];det=det+c[i]*w[i];}
    if(value(area2)<=0||value(det)<=0)throw std::runtime_error("Invalid directional crack orientation/deformation");
    T area=root(area2),q=uu+vv,sum=root(q+T(2)*area),height=det/area;
    T A=T(Real(4)/3*mu)+T(.5*K)*area2,B=T(Real(2)/3*mu)*sum;
    T free_height=(B+root(B*B+T(2*K)*A))/(T(2)*A);
    bool closed=value(height)<value(free_height);T elastic_height=closed?height:free_height;
    Real max_stretch=std::sqrt(.5*(value(q)+std::sqrt(std::max(Real(0),value(q)*value(q)-4*value(area2)))));
    Real min_stretch=value(area)/max_stretch;
    // The unchanged Biot solid is not globally rank-one convex. Do not hide a
    // nonconvex elastic shear branch behind a fake numerical regularization.
    if(max_stretch>=2*(min_stretch+value(elastic_height)))throw std::runtime_error("Directional crack relaxation left the stable finite-volume Biot shear branch");
    T anisotropy=((uu-vv)*(uu-vv)+T(4)*uv*uv)/(q+T(2)*area);
    T thickness=elastic_height-T(.5)*sum;
    T energy=T(mu)*(T(.5)*anisotropy+T(Real(2)/3)*thickness*thickness)+bulk(area*elastic_height,K);
    return {energy,height,elastic_height,closed};
}
}

namespace SPH::solid_dynamics
{
template<class Elastic>class DirectionalCrackEnergy
{
    const Elastic&elastic;
  public:
    struct Result{Real energy,tangent_bound,normal_height,elastic_height;Matd P;bool closed;};
    explicit DirectionalCrackEnergy(const Elastic&e):elastic(e){}
    Real energy(const Matd&F,const Vecd&normal)const
    {
        Vecd N=normal.normalized();std::array<Real,9>a{};for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j)a[3*i+j]=F(i,j);
        return directional_detail::evaluate(a,N,elastic.K,elastic.mu).energy;
    }
    Result evaluate(const Matd&F,const Vecd&normal)const
    {
        if(!normal.allFinite()||std::abs(normal.norm()-1)>1e-8)throw std::runtime_error("Directional crack requires a unit reference normal");
        using namespace directional_detail;std::array<Jet,9>a{};for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j)a[3*i+j]=Jet::variable(F(i,j),3*i+j);
        auto e=directional_detail::evaluate(a,normal,elastic.K,elastic.mu);Result r;r.energy=e.energy.v;r.normal_height=e.normal_height.v;r.elastic_height=e.elastic_height.v;r.closed=e.closed;r.tangent_bound=0;
        for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j)r.P(i,j)=e.energy.g[3*i+j];
        // Symmetric Hessian absolute block-row bound: conservative even in
        // finite rotations. It replaces the old inverse-F volumetric estimate.
        for(size_t i=0;i<9;++i){Real sum=0;for(size_t j=0;j<9;++j)sum+=std::abs(e.energy.h[9*i+j]);r.tangent_bound=std::max(r.tangent_bound,sum);}
        if(!std::isfinite(r.energy)||!r.P.allFinite()||!std::isfinite(r.tangent_bound))throw std::runtime_error("Nonfinite directional crack response");return r;
    }
    Real initiationEnergy(Real strength)const
    {
        Matd F=elastic.uniaxial(1);Real lo=1,hi=1+std::max(Real(.01),2*strength/elastic.E);
        auto traction=[&](Real a){auto e=elastic.evaluate(elastic.uniaxial(a));return e.positive_P(0,0)+e.negative_P(0,0);};
        for(int k=0;k<40&&traction(hi)<strength;++k)hi=1+2*(hi-1);
        if(traction(hi)<strength)throw std::runtime_error("Directional tensile calibration could not be bracketed");
        for(int k=0;k<80;++k){Real m=.5*(lo+hi);if(traction(m)>strength)hi=m;else lo=m;}
        F=elastic.uniaxial(.5*(lo+hi));auto e=elastic.evaluate(F);Real released=e.positive+e.negative-energy(F,Vecd::UnitX());
        if(released<=0)throw std::runtime_error("No directional tensile energy release at prescribed strength");return released;
    }
};
}
