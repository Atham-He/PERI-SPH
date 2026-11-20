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

class CudaCrackProjection final : public CrackProjectionBatch
{
  public:
    CudaCrackProjection()
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

        constexpr int threads = 128;
        const std::size_t blocks = (count + threads - 1) / threads;
        if (blocks > static_cast<std::size_t>(properties_.maxGridSize[0]))
            throw std::length_error("CUDA crack projection exceeds the device grid");
        input_.reserve(count);
        output_.reserve(count);
        check(cudaMemcpy(input_.get(), input.data(), count * sizeof(input[0]),
                         cudaMemcpyHostToDevice), "input transfer");
        projectBatch<<<static_cast<unsigned int>(blocks), threads>>>(input_.get(), output_.get(), count);
        check(cudaGetLastError(), "kernel launch");
        check(cudaMemcpy(output.data(), output_.get(), count * sizeof(output[0]),
                         cudaMemcpyDeviceToHost), "output transfer");
    }

  private:
    cudaDeviceProp properties_{};
    DeviceBuffer<CrackProjectionInput> input_;
    DeviceBuffer<CrackProjectionResult> output_;
};
}

std::unique_ptr<CrackProjectionBatch> makeCudaCrackProjection()
{
    return std::make_unique<CudaCrackProjection>();
}
}
