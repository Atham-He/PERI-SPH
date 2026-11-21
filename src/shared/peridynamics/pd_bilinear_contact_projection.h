#pragma once

#include "pd_crack_projection.h"
#include <cmath>
#include <limits>

#if defined(__CUDACC__)
#define SPH_CRACK_HOST_DEVICE __host__ __device__
#else
#define SPH_CRACK_HOST_DEVICE
#endif

namespace SPH::peridynamics
{
namespace crack_projection_detail
{
inline constexpr double face_tolerance = 1e-9;
inline constexpr double lower = -face_tolerance;
inline constexpr double upper = 1 + face_tolerance;
inline constexpr double epsilon = std::numeric_limits<double>::epsilon();

struct Vector3
{
    double x, y, z;
    SPH_CRACK_HOST_DEVICE Vector3 operator+(const Vector3 &r) const { return {x+r.x,y+r.y,z+r.z}; }
    SPH_CRACK_HOST_DEVICE Vector3 operator-(const Vector3 &r) const { return {x-r.x,y-r.y,z-r.z}; }
    SPH_CRACK_HOST_DEVICE Vector3 operator*(double s) const { return {s*x,s*y,s*z}; }
    SPH_CRACK_HOST_DEVICE double dot(const Vector3 &r) const { return x*r.x+y*r.y+z*r.z; }
    SPH_CRACK_HOST_DEVICE double norm() const { return std::sqrt(dot(*this)); }
    SPH_CRACK_HOST_DEVICE double absDot(const Vector3 &r) const
    {
        return std::abs(x*r.x)+std::abs(y*r.y)+std::abs(z*r.z);
    }
};

struct Face
{
    Vector3 edge_u, edge_v, warp, query;
    SPH_CRACK_HOST_DEVICE Vector3 tangentU(double t) const { return edge_u + warp*t; }
    SPH_CRACK_HOST_DEVICE Vector3 tangentV(double s) const { return edge_v + warp*s; }
    SPH_CRACK_HOST_DEVICE Vector3 residual(double s, double t) const
    {
        return query - (edge_u*s + edge_v*t + warp*(s*t));
    }
};

struct Quadratic
{
    double c[3];
};

SPH_CRACK_HOST_DEVICE inline Quadratic rulingLengthSquared(const Vector3 &edge, const Vector3 &warp)
{
    return {{edge.dot(edge), 2*edge.dot(warp), warp.dot(warp)}};
}

SPH_CRACK_HOST_DEVICE inline Quadratic rulingNumerator(const Vector3 &edge, const Vector3 &other,
                                                      const Vector3 &warp, const Vector3 &query)
{
    return {{query.dot(edge), query.dot(warp)-other.dot(edge), -other.dot(warp)}};
}

SPH_CRACK_HOST_DEVICE inline bool positiveOnUnitInterval(double c0, double c1, double c2, double margin)
{
    return c0 > margin && c0 + 0.5*c1 > margin && c0 + c1 + c2 > margin;
}

SPH_CRACK_HOST_DEVICE inline bool rulingExcludesFoot(const Vector3 &edge, const Vector3 &other,
                                                     const Vector3 &warp, const Vector3 &query)
{
    const auto a = rulingLengthSquared(edge, warp);
    const auto b = rulingNumerator(edge, other, warp, query);
    const double ruling_size = edge.norm() + warp.norm();
    const double margin = 16*face_tolerance*(query.norm() + other.norm() + ruling_size)*ruling_size;
    // These are sufficient sign tests on the quadratic Bernstein control
    // values of u*a(v)-b(v) at u=0 and u=1. The margin includes edge tolerance.
    return positiveOnUnitInterval(-b.c[0], -b.c[1], -b.c[2], margin) ||
           positiveOnUnitInterval(b.c[0]-a.c[0], b.c[1]-a.c[1], b.c[2]-a.c[2], margin);
}

SPH_CRACK_HOST_DEVICE inline bool inFace(double s, double t)
{
    return s >= lower && s <= upper && t >= lower && t <= upper;
}

SPH_CRACK_HOST_DEVICE inline bool localMinimum(const Face &face, double s, double t, double scale)
{
    const auto u = face.tangentU(t);
    const auto v = face.tangentV(s);
    const auto r = face.residual(s, t);
    const double a = u.dot(u), b = u.dot(v)-r.dot(face.warp), c = v.dot(v);
    const double tolerance = 512*epsilon*scale;
    return a > 0 && c > 0 && a*c >= b*b &&
           std::abs(u.dot(r)) <= tolerance*std::sqrt(a) &&
           std::abs(v.dot(r)) <= tolerance*std::sqrt(c);
}

SPH_CRACK_HOST_DEVICE inline void refineStationaryFoot(const Face &face, double scale, double &s, double &t)
{
    for (int iteration = 0; iteration < 40; ++iteration)
    {
        const auto u = face.tangentU(t);
        const auto v = face.tangentV(s);
        const auto r = face.residual(s, t);
        const double a = u.dot(u), c = v.dot(v), metric = std::fmax(a, c);
        if (!(metric > 0) || !std::isfinite(metric))
            return;
        const double gs = u.dot(r), gt = v.dot(r);
        if (std::abs(gs) <= 64*epsilon*scale*std::sqrt(a) &&
            std::abs(gt) <= 64*epsilon*scale*std::sqrt(c))
            return;

        double b = u.dot(v) - r.dot(face.warp);
        if (0.5*(a+c) - std::hypot(0.5*(a-c), b) <= 1e-13*metric)
            b = u.dot(v);
        if (a*c - b*b <= 1e-14*metric*metric)
            return;

        // Diagonally pivoted two-by-two LDLT, as in the original projection.
        double ds, dt;
        if (a >= c)
        {
            const double ratio = b/a;
            dt = (gt - ratio*gs)/(c - ratio*b);
            ds = gs/a - ratio*dt;
        }
        else
        {
            const double ratio = b/c;
            ds = (gs - ratio*gt)/(a - ratio*b);
            dt = gt/c - ratio*ds;
        }
        if (!std::isfinite(ds) || !std::isfinite(dt) || ds*ds + dt*dt < 4e-24)
            return;

        double alpha = 1;
        bool accepted = false;
        for (int line = 0; line < 25; ++line)
        {
            const double us = alpha*ds, vt = alpha*dt;
            if (inFace(s+us, t+vt))
            {
                const auto dx = u*us + v*vt + face.warp*(us*vt);
                const double norm_squared = dx.dot(dx);
                const double change = norm_squared - 2*r.dot(dx);
                if (change <= 64*epsilon*(norm_squared + 2*r.absDot(dx)))
                {
                    s += us;
                    t += vt;
                    accepted = true;
                    break;
                }
            }
            alpha *= 0.5;
        }
        if (!accepted)
            return;
    }
}

struct Polynomial
{
    double c[6];
    SPH_CRACK_HOST_DEVICE double value(int degree, double x) const
    {
        double result = c[degree];
        for (int k = degree-1; k >= 0; --k)
            result = result*x + c[k];
        return result;
    }
};

struct Roots
{
    // At most 61 trial roots for degree five, even if every endpoint that
    // rounds to zero is retained along with the adjacent sign-change roots.
    double x[61];
    int count;
};

SPH_CRACK_HOST_DEVICE inline double bisect(const Polynomial &p, int degree, double lo, double hi, double sign_lo)
{
    for (int step = 0; step < 54; ++step)
    {
        const double mid = 0.5*(lo+hi);
        const double value = p.value(degree, mid);
        if (value == 0 || mid == lo || mid == hi)
            return mid;
        if ((value < 0) == (sign_lo < 0))
            lo = mid;
        else
            hi = mid;
    }
    return 0.5*(lo+hi);
}

SPH_CRACK_HOST_DEVICE inline void isolateRoots(Polynomial polynomial, Roots &roots)
{
    int degree = 5;
    while (degree > 0 && polynomial.c[degree] == 0)
        --degree;
    roots.count = 0;
    if (degree == 0)
        return;

    Polynomial derivatives[6]{};
    derivatives[degree] = polynomial;
    for (int k = degree; k > 1; --k)
        for (int j = 1; j <= k; ++j)
            derivatives[k-1].c[j-1] = j*derivatives[k].c[j];

    const double linear = -derivatives[1].c[0]/derivatives[1].c[1];
    if (linear >= lower && linear <= upper)
        roots.x[roots.count++] = linear;
    for (int k = 2; k <= degree; ++k)
    {
        const auto &p = derivatives[k];
        double scale = 0;
        for (int j = 0; j <= k; ++j)
            scale += std::abs(p.c[j]);
        const double rounding = 16*epsilon*scale;
        Roots next;
        next.count = 0;
        double lo = lower, vlo = p.value(k, lo);
        for (int j = 0; j <= roots.count; ++j)
        {
            const double hi = j == roots.count ? upper : roots.x[j];
            const double vhi = p.value(k, hi);
            if (std::abs(vlo) <= rounding)
                next.x[next.count++] = lo;
            if (lo < hi && ((vlo < 0 && vhi > 0) || (vlo > 0 && vhi < 0)))
                next.x[next.count++] = bisect(p, k, lo, hi, vlo);
            lo = hi;
            vlo = vhi;
        }
        if (std::abs(vlo) <= rounding)
            next.x[next.count++] = upper;
        for (int j = 0; j < next.count; ++j)
            roots.x[j] = next.x[j];
        roots.count = next.count;
    }
}

SPH_CRACK_HOST_DEVICE inline CrackProjectionResult findRulingStationaryFeet(const Face &face)
{
    const auto a = rulingLengthSquared(face.edge_u, face.warp);
    const auto b = rulingNumerator(face.edge_u, face.edge_v, face.warp, face.query);
    const double dot_v = face.query.dot(face.edge_v), v_squared = face.edge_v.dot(face.edge_v);
    Polynomial derivative{};
    // With u=b(v)/a(v), the stationarity equation is the derivative of
    // the reduced squared distance multiplied by a(v)^2 / 2.
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            derivative.c[i+j] -= dot_v*a.c[i]*a.c[j];
            derivative.c[i+j+1] += v_squared*a.c[i]*a.c[j];
            for (int k = 1; k < 3; ++k)
                derivative.c[i+j+k-1] += k*b.c[i]*(0.5*b.c[j]*a.c[k] - a.c[j]*b.c[k]);
        }

    Roots roots;
    isolateRoots(derivative, roots);
    CrackProjectionResult best{0, 0, CrackProjectionResult::Outside};
    double distance = HUGE_VAL;
    const double scale = std::fmax(1.0, face.query.norm());
    for (int i = 0; i < roots.count; ++i)
    {
        double t = roots.x[i];
        const auto ruling = face.tangentU(t);
        const double length_squared = ruling.dot(ruling);
        if (length_squared <= 0)
            continue;
        double s = (face.query-face.edge_v*t).dot(ruling)/length_squared;
        if (!inFace(s, t))
            continue;
        refineStationaryFoot(face, scale, s, t);
        const auto r = face.residual(s, t);
        if (r.dot(r) < distance && localMinimum(face, s, t, scale))
        {
            best = {s, t, CrackProjectionResult::Accepted};
            distance = r.dot(r);
        }
    }
    return best;
}
}

/** Find an orthogonal contact foot on the finite bilinear face.
 * A boundary minimum caused by clamping an out-of-face point is not a
 * face contact in this law. Reference quadrature uses its own inversion.
 */
SPH_CRACK_HOST_DEVICE inline CrackProjectionResult projectCrackContact(const CrackProjectionInput &input)
{
    using namespace crack_projection_detail;
    for (int i = 0; i < 5; ++i)
        for (int d = 0; d < 3; ++d)
            if (!std::isfinite(input.coordinates[i][d]))
                return {0, 0, CrackProjectionResult::InvalidGeometry};
    const auto &x = input.coordinates;
    const Vector3 origin{x[0][0], x[0][1], x[0][2]};
    const Vector3 u = Vector3{x[1][0], x[1][1], x[1][2]}-origin;
    const Vector3 diagonal = Vector3{x[2][0], x[2][1], x[2][2]}-origin;
    const Vector3 v = Vector3{x[3][0], x[3][1], x[3][2]}-origin;
    const Face face{u, v, (diagonal-v)-u, Vector3{x[4][0], x[4][1], x[4][2]}-origin};
    const double scale = std::fmax(u.norm(), std::fmax(v.norm(), face.warp.norm()));
    if (!(scale > 0) || !std::isfinite(scale) || !std::isfinite(face.query.norm()))
        return {0, 0, CrackProjectionResult::InvalidGeometry};

    if (rulingExcludesFoot(u, v, face.warp, face.query) || rulingExcludesFoot(v, u, face.warp, face.query))
        return {0, 0, CrackProjectionResult::Outside};

    double s = 0.5, t = 0.5;
    const double scale_with_query = std::fmax(scale, face.query.norm());
    refineStationaryFoot(face, scale_with_query, s, t);
    if (localMinimum(face, s, t, scale_with_query))
        return {s, t, CrackProjectionResult::Accepted};

    const double inv = 1/scale;
    return findRulingStationaryFeet({u*inv, v*inv, face.warp*inv, face.query*inv});
}
}

#undef SPH_CRACK_HOST_DEVICE
