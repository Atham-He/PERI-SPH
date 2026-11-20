#include "pd_crack_projection.h"
#include "pd_cuda_buffer.cuh"

#include <cuda_runtime.h>
#include <cmath>

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;

// The double-precision solve follows BilinearCrackGeometry::closest,
// including its stationarity, damping and rank checks.
struct Vector3 {
    double x, y, z;
    __device__ Vector3 operator+(const Vector3& r) const {return {x+r.x,y+r.y,z+r.z};}
    __device__ Vector3 operator-(const Vector3& r) const {return {x-r.x,y-r.y,z-r.z};}
    __device__ Vector3 operator*(double s) const {return {s*x,s*y,s*z};}
    __device__ double dot(const Vector3& r) const {return x*r.x+y*r.y+z*r.z;}
    __device__ double norm() const {return sqrt(dot(*this));}
    __device__ double absDot(const Vector3& r) const {return fabs(x*r.x)+fabs(y*r.y)+fabs(z*r.z);}
};

__device__ CrackProjectionResult projectOne(const CrackProjectionInput &data) {
    constexpr double eps=2.2204460492503131e-16;
    const Vector3 origin{data.coordinates[0][0],data.coordinates[0][1],data.coordinates[0][2]};
    const Vector3 e=Vector3{data.coordinates[1][0],data.coordinates[1][1],data.coordinates[1][2]}-origin;
    const Vector3 diagonal=Vector3{data.coordinates[2][0],data.coordinates[2][1],data.coordinates[2][2]}-origin;
    const Vector3 f=Vector3{data.coordinates[3][0],data.coordinates[3][1],data.coordinates[3][2]}-origin;
    const Vector3 q=Vector3{data.coordinates[4][0],data.coordinates[4][1],data.coordinates[4][2]}-origin;
    const Vector3 warp=(diagonal-f)-e;
    const double coordinate_scale=fmax(fmax(e.norm(),f.norm()),fmax(q.norm(),fmax(diagonal.norm(),1.e-30)));
    double s=.5,t=.5;
    for (int iteration=0;iteration<40;++iteration) {
        const Vector3 u=e+warp*t,v=f+warp*s;
        const Vector3 residual=q-(e*s+f*t+warp*(s*t));
        const double a=u.dot(u),c=v.dot(v),scale=fmax(a,c);
        if (!(scale>0) || !isfinite(scale)) return {s,t,1};
        const double gs=u.dot(residual),gt=v.dot(residual);
        const double rounding=64*eps*coordinate_scale;
        if (fabs(gs)<=rounding*sqrt(a) && fabs(gt)<=rounding*sqrt(c)) return {s,t,0};
        double b=u.dot(v)-residual.dot(warp);
        const double mean=.5*(a+c), radius=hypot(.5*(a-c),b);
        if (mean-radius<=1.e-13*scale) b=u.dot(v);
        if (a*c-b*b<=1.e-14*scale*scale) return {s,t,2};
        double ds,dt;
        if (a>=c) {
            const double l=b/a;
            dt=(gt-l*gs)/(c-l*b);
            ds=gs/a-l*dt;
        } else {
            const double l=b/c;
            ds=(gs-l*gt)/(a-l*b);
            dt=gt/c-l*ds;
        }
        if (!isfinite(ds) || !isfinite(dt)) return {s,t,3};
        if (sqrt(ds*ds+dt*dt)<2.e-12) return {s,t,0};
        double alpha=1;
        bool accepted=false;
        for (int line=0;line<25;++line) {
            const double u_step=alpha*ds,v_step=alpha*dt;
            const Vector3 dx=u*u_step+v*v_step+warp*(u_step*v_step);
            const double norm2=dx.dot(dx),change=norm2-2*residual.dot(dx);
            const double tolerance=64*eps*(norm2+2*residual.absDot(dx));
            if (change<=tolerance) {s+=u_step;t+=v_step;accepted=true;break;}
            alpha*=.5;
        }
        if (!accepted) return {s,t,sqrt(gs*gs+gt*gt)<1.e-12*scale?0:4};
    }
    return {s,t,5};
}

__global__ void projectBatch(const CrackProjectionInput *in, CrackProjectionResult *out,std::size_t n) {
    const std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) out[i] = projectOne(in[i]);
}

__device__ Vector3 loadVector(const double (&x)[3])
{
    return {x[0], x[1], x[2]};
}

__device__ CrackLocalContactResult localContact(const CrackLocalContactInput &input)
{
    CrackLocalContactResult response{};
    const Vector3 u = loadVector(input.tangent_u);
    const Vector3 v = loadVector(input.tangent_v);
    const Vector3 normal = loadVector(input.normal);
    const double penetration = fmin(0.0, input.gap);
    const Vector3 residual = normal * penetration;
    const double a = u.dot(u);
    const double b = u.dot(v) - residual.dot(loadVector(input.mixed_tangent));
    const double c = v.dot(v);
    const double largest = .5 * (a + c) + hypot(.5 * (a - c), b);
    const double determinant = a * c - b * b;

    // det / largest is the smaller eigenvalue of the symmetric tangent.
    // The host eigensolver handles borderline geometry and reports failures.
    if (!isfinite(determinant) || !isfinite(largest) || largest <= 0.0 ||
        determinant <= 1.e-13 * largest * largest)
        return response;

    const double kinv00 = c / determinant;
    const double kinv01 = -b / determinant;
    const double kinv11 = a / determinant;
    double coefficient[5] = {1.0, -input.shape[0], -input.shape[1], -input.shape[2], -input.shape[3]};
    double coupling_u[15], coupling_v[15];
    for (int i = 0; i < 5; ++i)
    {
        const Vector3 pu = i == 0 ? u * (-coefficient[i]) :
            u * (-coefficient[i]) - residual * input.derivative_u[i - 1];
        const Vector3 pv = i == 0 ? v * (-coefficient[i]) :
            v * (-coefficient[i]) - residual * input.derivative_v[i - 1];
        coupling_u[3 * i] = pu.x;
        coupling_u[3 * i + 1] = pu.y;
        coupling_u[3 * i + 2] = pu.z;
        coupling_v[3 * i] = pv.x;
        coupling_v[3 * i + 1] = pv.y;
        coupling_v[3 * i + 2] = pv.z;
    }

    for (int row = 0; row < 15; ++row)
    {
        const double inverse_u = -coupling_u[row] * kinv00 - coupling_v[row] * kinv01;
        const double inverse_v = -coupling_u[row] * kinv01 - coupling_v[row] * kinv11;
        double bound = 0.0;
        for (int column = 0; column < 15; ++column)
        {
            // Condense out the motion of the two closest-point coordinates.
            double h = inverse_u * coupling_u[column] + inverse_v * coupling_v[column];
            if (row % 3 == column % 3)
                h += coefficient[row / 3] * coefficient[column / 3];
            bound += fabs(input.stiffness * h) * input.mass_weight[column / 3];
        }
        if (!isfinite(bound))
            return response;
        response.row_bound[row / 3][row % 3] = bound;
    }
    if (input.gap < 0.0)
    {
        const Vector3 force = normal * (-input.stiffness * input.gap);
        response.force_on_slave[0] = force.x;
        response.force_on_slave[1] = force.y;
        response.force_on_slave[2] = force.z;
    }
    response.valid = true;
    return response;
}

__global__ void localContactBatch(const CrackLocalContactInput *input, CrackLocalContactResult *output,
                                 std::size_t count)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < count)
        output[i] = localContact(input[i]);
}

class CudaCrackProjection final : public CrackProjectionBatch
{
  public:
    explicit CudaCrackProjection(bool assemble_local_contact) : assemble_local_contact_(assemble_local_contact)
    {
        int device = 0;
        check(cudaGetDevice(&device), "device selection");
        check(cudaGetDeviceProperties(&properties_, device), "device properties");
    }

    void project(const std::vector<CrackProjectionInput> &input,
                 std::vector<CrackProjectionResult> &output) override
    {
        const std::size_t count = input.size();
        output.resize(count);
        if (count == 0)
            return;

        const auto blocks = grid(count);
        input_.reserve(count);
        output_.reserve(count);
        check(cudaMemcpy(input_.get(), input.data(), count * sizeof(input[0]),
                         cudaMemcpyHostToDevice), "input transfer");
        projectBatch<<<blocks, threads>>>(input_.get(), output_.get(), count);
        check(cudaGetLastError(), "kernel launch");
        check(cudaMemcpy(output.data(), output_.get(), count * sizeof(output[0]),
                         cudaMemcpyDeviceToHost), "output transfer");
    }

    bool assemblesLocalContact() const noexcept override { return assemble_local_contact_; }

    void assembleLocalContact(const std::vector<CrackLocalContactInput> &input,
                              std::vector<CrackLocalContactResult> &output) override
    {
        if (!assemble_local_contact_)
            throw std::logic_error("Local contact assembly is disabled for this projection backend");
        const auto count = input.size();
        output.resize(count);
        if (count == 0)
            return;
        const auto blocks = grid(count);
        contact_input_.reserve(count);
        contact_output_.reserve(count);
        check(cudaMemcpy(contact_input_.get(), input.data(), count * sizeof(input[0]),
                         cudaMemcpyHostToDevice), "crack tangent input transfer");
        localContactBatch<<<blocks, threads>>>(contact_input_.get(), contact_output_.get(), count);
        check(cudaGetLastError(), "crack tangent kernel launch");
        check(cudaMemcpy(output.data(), contact_output_.get(), count * sizeof(output[0]),
                         cudaMemcpyDeviceToHost), "crack tangent output transfer");
    }

  private:
    unsigned int grid(std::size_t count) const
    {
        const auto blocks = (count + threads - 1) / threads;
        if (blocks > static_cast<std::size_t>(properties_.maxGridSize[0]))
            throw std::length_error("CUDA crack batch exceeds the device grid");
        return static_cast<unsigned int>(blocks);
    }

    static constexpr unsigned int threads = 128;
    bool assemble_local_contact_;
    cudaDeviceProp properties_{};
    DeviceBuffer<CrackProjectionInput> input_;
    DeviceBuffer<CrackProjectionResult> output_;
    DeviceBuffer<CrackLocalContactInput> contact_input_;
    DeviceBuffer<CrackLocalContactResult> contact_output_;
};
}

std::unique_ptr<CrackProjectionBatch> makeCudaCrackProjection(bool assemble_local_contact)
{
    return std::make_unique<CudaCrackProjection>(assemble_local_contact);
}
}
