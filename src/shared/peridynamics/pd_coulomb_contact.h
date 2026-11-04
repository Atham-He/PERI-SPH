/* Licensed under the Apache License, Version 2.0. */
#ifndef PD_COULOMB_CONTACT_H
#define PD_COULOMB_CONTACT_H

#include "pd_many_contact.h"
#include "pd_coulomb_disk.h"

namespace SPH::peridynamics
{
using ContactSpatialVector = Eigen::Matrix<Real,3,1>;

/** All three gather rows are built from ONE signed partition and frame.
 * Impulses are point impulses, not quadrature pressure densities.
 * The wall tangent velocity is the caller's declared actual-interval cursor.
 */
struct CoulombContactConstraint
{
    StdVec<Real> weights;
    ContactSpatialVector normal = ContactSpatialVector::UnitX();
    std::array<ContactSpatialVector,2> tangent{
        ContactSpatialVector::UnitY(), ContactSpatialVector::UnitZ()};
    Real old_offset=0, new_offset=0, mu=0;
    ContactTangentVector wall_tangent_velocity=ContactTangentVector::Zero();
};

namespace coulomb_contact_detail
{
using Rational=ContactExactScalar;
using Integer=exact_contact_detail::Integer;
using Wide=coulomb_disk_detail::Wide;

inline Rational exactReal(Real x)
{
    const auto binary=exact_contact_detail::decode(x);
    Rational result(binary.mantissa);
    if (binary.exponent>=0) result*=Integer(1)<<binary.exponent;
    else result/=Integer(1)<<(-binary.exponent);
    return result;
}

/** Accumulate binary64 products as integers before forming ONE rational.
 * This avoids both premature binary64 gather rounding and repeated rational
 * normalization per product. It does not alter any supplied coefficient.
 */
inline Rational exactDot(const StdVec<Real> &a,const StdVec<Real> &b)
{
    if (a.size()!=b.size() || a.size()>32768)
        throw std::invalid_argument("Coulomb contact: exact dot dimensions");
    Integer sum=0;
    for (size_t i=0;i<a.size();++i)
    {
        const auto x=exact_contact_detail::decode(a[i]),y=exact_contact_detail::decode(b[i]);
        if (x.mantissa!=0 && y.mantissa!=0)
            sum+=(x.mantissa*y.mantissa)<<(x.exponent+y.exponent+2148);
    }
    return Rational(sum)/Rational(Integer(1)<<2148);
}

inline Rational absolute(const Rational &x) { return x<0 ? -x : x; }
inline Rational positive(const Rational &x) { return x>0 ? x : Rational(0); }

inline Real rounded(const Rational &value,Real &bound,CoulombDiskValue &certificate,
    const Wide &additional_error=Wide(0))
{
    const Wide wide=contactExactToWide(value);
    return coulomb_disk_detail::convert(wide,
        coulomb_disk_detail::arithmeticEnvelope(coulomb_disk_detail::absolute(wide))+additional_error,
        bound,certificate);
}

inline Real upper(const Rational &value,CoulombDiskValue &certificate)
{
    if (value<0) throw std::logic_error("Coulomb contact: negative upper metric");
    const Wide wide=contactExactToWide(value);
    return coulomb_disk_detail::outwardPositive(wide+
        coulomb_disk_detail::arithmeticEnvelope(wide),certificate);
}

inline bool frameValid(const CoulombContactConstraint &row)
{
    const std::array<ContactSpatialVector,3> basis{row.normal,row.tangent[0],row.tangent[1]};
    for (const auto &v:basis) if (!v.allFinite()) return false;
    const Rational tolerance=exactReal(1e-12);
    for (size_t i=0;i<3;++i) for (size_t j=0;j<3;++j)
    {
        StdVec<Real> a{basis[i][0],basis[i][1],basis[i][2]},b{basis[j][0],basis[j][1],basis[j][2]};
        if (absolute(exactDot(a,b)-Rational(i==j))>tolerance) return false;
    }
    return true;
}
} // namespace coulomb_contact_detail

struct CoulombContactValues
{
    StdVec<Real> residual,residual_roundoff_bound,supplied_momentum;
    StdVec<Real> normalized_gaps,gap_roundoff_bound;
    StdVec<ContactTangentVector> normalized_relative_velocity;
    StdVec<ContactTangentVector> relative_velocity_roundoff_bound;
    StdVec<CoulombDiskValue> disk;
    DenseContactMatrix jacobian,supplied_jacobian;
    Real merit=0,squared_norm=0,KKT=0,maximum_jacobian_roundoff_bound=0;
    bool precision_budget_ok=true;
};

/** Immutable original coefficients, geometry bounds and normalization.
 * No reference/history is refreshed by nonlinear trials. The caller supplies
 * the original material momentum/Jacobian and is the only possible committer.
 */
class CoulombContactSystem
{
  public:
    const size_t dofs,contacts_count,unknowns;
    const Real dt,speed,impulse_scale;
    const StdVec<CoulombContactConstraint> contacts;
    const DenseContactMatrix &coefficients() const { return gather; }

    CoulombContactSystem(const StdVec<Real> &old_position,Real actual_dt,Real mass_scale,
        const StdVec<CoulombContactConstraint> &rows,const ManyContactSettings &settings={})
        :dofs(old_position.size()),contacts_count(rows.size()),unknowns(dofs+3*contacts_count),
         dt(actual_dt),speed(settings.velocity_scale),impulse_scale(mass_scale*speed),contacts(rows)
    {
        using namespace coulomb_contact_detail;
        if (!dofs || dofs%3 || dofs>128 || contacts_count>256 ||
            !std::isfinite(dt) || dt<=0 || !std::isfinite(mass_scale) || mass_scale<=0 ||
            !std::isfinite(speed) || speed<=0 || !std::isfinite(impulse_scale) || impulse_scale<=0 ||
            !std::isfinite(settings.relative_tolerance) || settings.relative_tolerance<=0 ||
            settings.relative_tolerance>1e-12 || !settings.max_newton_steps || settings.max_newton_steps>48 ||
            !settings.max_line_search || settings.max_line_search>32)
            throw std::invalid_argument("Coulomb contact: bounded dimensions/settings");
        for (Real x:old_position) if (!std::isfinite(x))
            throw std::invalid_argument("Coulomb contact: nonfinite accepted position");
        gather.assign(3*contacts_count,StdVec<Real>(dofs));
        normal_bound_.resize(contacts_count);
        const Rational exact_dt=exactReal(dt),exact_speed=exactReal(speed);
        for (size_t s=0;s<contacts_count;++s)
        {
            const auto &row=contacts[s];
            if (row.weights.size()!=dofs/3 || !frameValid(row) ||
                !std::isfinite(row.old_offset) || !std::isfinite(row.new_offset) ||
                !std::isfinite(row.mu) || row.mu<0 || !row.wall_tangent_velocity.allFinite())
                throw std::invalid_argument("Coulomb contact: invalid row/frame/coefficient");
            Rational partition=0;
            for (Real w:row.weights) partition+=exactReal(w);
            if (absolute(partition-1)>exactReal(1e-12))
                throw std::invalid_argument("Coulomb contact: signed partition mismatch");
            const std::array<ContactSpatialVector,3> basis{row.normal,row.tangent[0],row.tangent[1]};
            for (size_t a=0;a<3;++a) for (size_t node=0;node<row.weights.size();++node)
                for (size_t component=0;component<3;++component)
                {
                    const Real coefficient=row.weights[node]*basis[a][component];
                    if (!std::isfinite(coefficient))
                        throw std::invalid_argument("Coulomb contact: nonfinite constructed coefficient");
                    gather[3*s+a][3*node+component]=coefficient;
                }
            const Rational old_gather=exactDot(gather[3*s],old_position);
            if ((old_gather-exactReal(row.old_offset))/(exact_dt*exact_speed)<-exactReal(settings.relative_tolerance))
                throw std::invalid_argument("Coulomb contact: accepted initial penetration");
            normal_bound_[s]=(exactReal(row.new_offset)-old_gather)/(exact_dt*exact_speed);
        }
    }

    CoulombContactValues evaluate(const StdVec<Real> &velocity,const StdVec<Real> &impulses,
        const std::function<NonlinearContactEquation(const StdVec<Real> &)> &material,bool tangent=true) const
    {
        using namespace coulomb_contact_detail;
        if (velocity.size()!=dofs || impulses.size()!=3*contacts_count)
            throw std::invalid_argument("Coulomb contact: iterate dimensions");
        for (Real x:velocity) if (!std::isfinite(x)) throw std::invalid_argument("Coulomb contact: nonfinite velocity");
        for (Real x:impulses) if (!std::isfinite(x)) throw std::invalid_argument("Coulomb contact: nonfinite impulse");
        const auto equation=material(velocity);
        if (equation.momentum.size()!=dofs || equation.jacobian.size()!=dofs)
            throw std::runtime_error("Coulomb contact: material callback dimensions");
        CoulombContactValues out;out.residual.assign(unknowns,0);out.residual_roundoff_bound.assign(unknowns,0);
        out.normalized_gaps.resize(contacts_count);out.gap_roundoff_bound.resize(contacts_count);
        out.normalized_relative_velocity.resize(contacts_count);out.relative_velocity_roundoff_bound.resize(contacts_count);
        out.disk.resize(contacts_count);out.supplied_momentum=equation.momentum;out.supplied_jacobian=equation.jacobian;
        if (tangent) out.jacobian.assign(unknowns,StdVec<Real>(unknowns,0));
        CoulombDiskValue arithmetic;const Rational I=exactReal(impulse_scale),V=exactReal(speed);
        auto jacobianValue=[&](size_t i,size_t j,const Rational &value,const Wide &additional_error=Wide(0))
        {
            if (!tangent) return;
            Real bound=0;out.jacobian[i][j]=rounded(value,bound,arithmetic,additional_error);
            out.maximum_jacobian_roundoff_bound=std::max(out.maximum_jacobian_roundoff_bound,bound);
        };
        for (size_t i=0;i<dofs;++i)
        {
            if (equation.jacobian[i].size()!=dofs || !std::isfinite(equation.momentum[i]))
                throw std::runtime_error("Coulomb contact: material callback row");
            StdVec<Real> coefficients(impulses.size());
            for (size_t s=0;s<coefficients.size();++s) coefficients[s]=gather[s][i];
            const auto balance=exactContactBalance(equation.momentum[i],coefficients,impulses,impulse_scale);
            out.residual[i]=balance.normalized_value;out.residual_roundoff_bound[i]=balance.normalized_roundoff_bound;
            out.precision_budget_ok &= balance.within_budget;
            for (size_t j=0;j<dofs;++j)
            {
                if (!std::isfinite(equation.jacobian[i][j]))
                    throw std::runtime_error("Coulomb contact: nonfinite original material tangent");
                jacobianValue(i,j,exactReal(equation.jacobian[i][j])*V/I);
            }
            if (tangent) for (size_t s=0;s<coefficients.size();++s) out.jacobian[i][dofs+s]=-coefficients[s];
        }
        Rational normal_products=0;
        Wide dissipation_errors=0;
        for (size_t s=0;s<contacts_count;++s)
        {
            const size_t normal_index=dofs+3*s;
            const Rational l=exactReal(impulses[3*s])/I;
            const Rational gap=exactDot(gather[3*s],velocity)/V-normal_bound_[s];
            out.normalized_gaps[s]=rounded(gap,out.gap_roundoff_bound[s],arithmetic);
            normal_products+=absolute(l*gap);
            out.KKT=std::max(out.KKT,upper(positive(-l),arithmetic));
            out.KKT=std::max(out.KKT,upper(positive(-gap),arithmetic));
            const Wide a=contactExactToWide(l),b=contactExactToWide(gap);
            const Wide radius=boost::multiprecision::sqrt(a*a+b*b);
            const Wide phi=radius-a-b;
            const Wide envelope=coulomb_disk_detail::arithmeticEnvelope(coulomb_disk_detail::absolute(a)+coulomb_disk_detail::absolute(b));
            out.residual[normal_index]=coulomb_disk_detail::convert(phi,envelope,
                out.residual_roundoff_bound[normal_index],arithmetic);
            const Wide da=radius==0 ? boost::multiprecision::sqrt(Wide(.5))-1 : a/radius-1;
            const Wide db=radius==0 ? da : b/radius-1;
            if (tangent)
            {
                Real derivative_bound=0;
                out.jacobian[normal_index][normal_index]=coulomb_disk_detail::convert(da,
                    coulomb_disk_detail::arithmeticEnvelope(2),derivative_bound,arithmetic);
                out.maximum_jacobian_roundoff_bound=std::max(out.maximum_jacobian_roundoff_bound,derivative_bound);
                for (size_t j=0;j<dofs;++j)
                {
                    const Wide value=db*Wide(gather[3*s][j]);
                    out.jacobian[normal_index][j]=coulomb_disk_detail::convert(value,
                        coulomb_disk_detail::arithmeticEnvelope(coulomb_disk_detail::absolute(value))+
                        coulomb_disk_detail::absolute(Wide(gather[3*s][j]))*coulomb_disk_detail::arithmeticEnvelope(2),
                        derivative_bound,arithmetic);
                    out.maximum_jacobian_roundoff_bound=std::max(out.maximum_jacobian_roundoff_bound,derivative_bound);
                }
            }
            std::array<Rational,2> t,w;
            for (size_t a=0;a<2;++a)
            {
                t[a]=exactReal(impulses[3*s+1+a])/I;
                w[a]=(exactDot(gather[3*s+1+a],velocity)-exactReal(contacts[s].wall_tangent_velocity[a]))/V;
                out.normalized_relative_velocity[s][a]=rounded(w[a],out.relative_velocity_roundoff_bound[s][a],arithmetic);
            }
            const auto disk=coulombDiskNaturalMapExact(l,t,w,exactReal(contacts[s].mu));out.disk[s]=disk;
            out.precision_budget_ok &= disk.within_budget;
            out.KKT=std::max(out.KKT,disk.cone_excess_upper);
            dissipation_errors+=Wide(disk.maximum_dissipation_error_upper);
            for (size_t a=0;a<2;++a)
            {
                const size_t row=normal_index+1+a;
                out.residual[row]=disk.residual[a];out.residual_roundoff_bound[row]=disk.residual_roundoff_bound[a];
                if (!tangent) continue;
                out.jacobian[row][normal_index]=disk.normal_derivative[a];
                out.maximum_jacobian_roundoff_bound=std::max(out.maximum_jacobian_roundoff_bound,disk.maximum_derivative_roundoff_bound);
                for (size_t b=0;b<2;++b) out.jacobian[row][normal_index+1+b]=disk.impulse_derivative(a,b);
                for (size_t j=0;j<dofs;++j)
                {
                    const Rational value=exactReal(disk.velocity_derivative(a,0))*exactReal(gather[3*s+1][j])+
                        exactReal(disk.velocity_derivative(a,1))*exactReal(gather[3*s+2][j]);
                    const Wide input_error=Wide(disk.maximum_derivative_roundoff_bound)*
                        (coulomb_disk_detail::absolute(Wide(gather[3*s+1][j]))+coulomb_disk_detail::absolute(Wide(gather[3*s+2][j])));
                    jacobianValue(row,j,value,input_error);
                }
            }
        }
        out.KKT=std::max(out.KKT,upper(normal_products,arithmetic));
        out.KKT=std::max(out.KKT,coulomb_disk_detail::outwardPositive(dissipation_errors+
            coulomb_disk_detail::arithmeticEnvelope(dissipation_errors),arithmetic));
        out.precision_budget_ok &= arithmetic.within_budget;
        if (!out.precision_budget_ok) return out;
        for (size_t i=0;i<unknowns;++i)
        {
            out.merit=std::max(out.merit,std::abs(out.residual[i])+out.residual_roundoff_bound[i]);
            out.squared_norm+=out.residual[i]*out.residual[i];
        }
        if (!std::isfinite(out.merit) || !std::isfinite(out.squared_norm) || !std::isfinite(out.KKT))
            throw std::runtime_error("Coulomb contact: nonfinite original merit/KKT");
        return out;
    }

  private:
    DenseContactMatrix gather; // normal,tangent0,tangent1 for each contact
    StdVec<coulomb_contact_detail::Rational> normal_bound_;
};

struct CoulombContactTrace
{
    size_t iteration=0,linear_rank=0,domain_rejections=0;
    Real alpha=0,linear_residual=0;
    StdVec<Real> velocity,impulses,normalized_direction;
    CoulombContactValues original;
};

struct CoulombContactResult
{
    bool converged=false;
    std::string status;
    StdVec<Real> velocity,impulses; // Empty unless every original equation/KKT passes.
    StdVec<CoulombContactTrace> trace;
    Real normalized_KKT=0;
    struct PrecisionWitness
    {
        bool present=false,during_line_search=false;
        size_t iteration=0;
        Real alpha=0;
        StdVec<Real> velocity,impulses;
        CoulombContactValues original;
    } precision_rejection;
    struct EvaluationFailure
    {
        bool present=false,during_line_search=false;
        size_t iteration=0;
        Real alpha=0;
        StdVec<Real> velocity,impulses;
        std::string reason;
    } evaluation_failure;
};

/** One nonsymmetric all-contact semismooth Newton method. No equation shifts,
 * row selection, gradient fallback, reference refresh or material commit.
 */
inline CoulombContactResult solveCoulombContact(const CoulombContactSystem &system,
    const std::function<NonlinearContactEquation(const StdVec<Real> &)> &material,
    const std::function<bool(const StdVec<Real> &)> &in_domain,const ManyContactSettings &settings={})
{
    if (settings.velocity_scale!=system.speed || !std::isfinite(settings.relative_tolerance) ||
        settings.relative_tolerance<=0 || settings.relative_tolerance>1e-12 ||
        !settings.max_newton_steps || settings.max_newton_steps>48 || !settings.max_line_search || settings.max_line_search>32)
        throw std::invalid_argument("Coulomb contact: solver settings mismatch");
    CoulombContactResult result;
    StdVec<Real> velocity(system.dofs,0),impulses(3*system.contacts_count,0);
    if (!in_domain(velocity)) throw std::invalid_argument("Coulomb contact: invalid initial iterate domain");
    auto reject=[&](CoulombContactValues original,const StdVec<Real> &v,const StdVec<Real> &p,
        size_t iteration,Real alpha,bool line_search)
    {
        result.status="original_contact_arithmetic_precision_budget";
        auto &w=result.precision_rejection;w.present=true;w.during_line_search=line_search;
        w.iteration=iteration;w.alpha=alpha;w.velocity=v;w.impulses=p;w.original=std::move(original);
    };
    auto failed_evaluation=[&](const std::exception &error,const StdVec<Real> &v,const StdVec<Real> &p,
        size_t iteration,Real alpha,bool line_search)
    {
        result.status="original_material_or_contact_evaluation_failed";
        auto &w=result.evaluation_failure;w.present=true;w.during_line_search=line_search;
        w.iteration=iteration;w.alpha=alpha;w.velocity=v;w.impulses=p;w.reason=error.what();
    };
    for (size_t iteration=0;iteration<=settings.max_newton_steps;++iteration)
    {
        CoulombContactValues current;
        try { current=system.evaluate(velocity,impulses,material,true); }
        catch(const std::exception &error)
        { failed_evaluation(error,velocity,impulses,iteration,0,false);return result; }
        if (!current.precision_budget_ok) { reject(std::move(current),velocity,impulses,iteration,0,false);return result; }
        CoulombContactTrace trace;trace.iteration=iteration;trace.velocity=velocity;trace.impulses=impulses;
        trace.original=std::move(current);
        if (std::max(trace.original.merit,trace.original.KKT)<=settings.relative_tolerance)
        {
            result.converged=true;result.status="all_normal_and_two_tangent_original_KKT_accepted";
            result.velocity=velocity;result.impulses=impulses;result.normalized_KKT=std::max(trace.original.merit,trace.original.KKT);
            result.trace.push_back(std::move(trace));return result;
        }
        if (iteration==settings.max_newton_steps)
        { result.trace.push_back(std::move(trace));result.status="Newton_iteration_limit";return result; }
        StdVec<Real> rhs(system.unknowns);
        for (size_t i=0;i<system.unknowns;++i) rhs[i]=-trace.original.residual[i];
        ConsistentLinearSolution direction;
        try { direction=solveBoundedConsistentContactSystem(trace.original.jacobian,rhs,1024); }
        catch (const std::runtime_error &)
        { result.trace.push_back(std::move(trace));result.status="original_coupled_linear_system_rejected";return result; }
        trace.linear_rank=direction.rank;trace.linear_residual=direction.relative_residual;trace.normalized_direction=direction.value;
        Real alpha=1;bool accepted=false;
        for (size_t search=0;search<settings.max_line_search;++search,alpha*=.5)
        {
            auto v=velocity,p=impulses;
            for (size_t i=0;i<system.dofs;++i) v[i]+=alpha*system.speed*direction.value[i];
            for (size_t s=0;s<p.size();++s) p[s]+=alpha*system.impulse_scale*direction.value[system.dofs+s];
            if (!in_domain(v)) { ++trace.domain_rejections;continue; }
            CoulombContactValues trial;
            try { trial=system.evaluate(v,p,material,false); }
            catch(const std::exception &error)
            {
                result.trace.push_back(std::move(trace));
                failed_evaluation(error,v,p,iteration,alpha,true);return result;
            }
            if (!trial.precision_budget_ok)
            {
                result.trace.push_back(std::move(trace));reject(std::move(trial),v,p,iteration,alpha,true);return result;
            }
            if (trial.squared_norm<=(1-1e-4*alpha)*trace.original.squared_norm)
            { trace.alpha=alpha;velocity=std::move(v);impulses=std::move(p);accepted=true;break; }
        }
        result.trace.push_back(std::move(trace));
        if (!accepted) { result.status="line_search_limit";return result; }
    }
    throw std::logic_error("Coulomb contact: unreachable");
}
} // namespace SPH::peridynamics
#endif
