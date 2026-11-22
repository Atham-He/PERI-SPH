#include "pd_gel_kernel.h"
#include "pd_cuda_buffer.cuh"

namespace SPH::peridynamics
{
namespace
{
using cuda_detail::check;
using cuda_detail::DeviceBuffer;
using cuda_detail::download;
using cuda_detail::upload;
// The device uses fixed-size arithmetic for the same constitutive kernel.
struct Matrix
{
    using Scalar = double;
    static constexpr int RowsAtCompileTime = 3, ColsAtCompileTime = 3;
    double data[3][3];

    __host__ __device__ double operator()(int a, int b) const { return data[a][b]; }
    __host__ __device__ double &operator()(int a, int b) { return data[a][b]; }

    __host__ __device__ static Matrix Identity()
    {
        Matrix out{};
        for (int i = 0; i < 3; ++i)
            out(i,i) = 1;
        return out;
    }

    __host__ __device__ Matrix operator+(const Matrix &other) const
    {
        Matrix out;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                out(a,b) = data[a][b] + other(a,b);
        return out;
    }

    __host__ __device__ Matrix operator-(const Matrix &other) const
    {
        Matrix out;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                out(a,b) = data[a][b] - other(a,b);
        return out;
    }

    __host__ __device__ Matrix scale(double factor) const
    {
        Matrix out;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                out(a,b) = factor * data[a][b];
        return out;
    }

    __host__ __device__ friend Matrix operator*(double factor, const Matrix &tensor) { return tensor.scale(factor); }
    __host__ __device__ Matrix operator/(double divisor) const { return scale(1.0 / divisor); }

    __host__ __device__ Matrix operator*(const Matrix &other) const
    {
        Matrix out{};
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                for (int k = 0; k < 3; ++k)
                    out(a,b) += data[a][k] * other(k,b);
        return out;
    }

    __host__ __device__ Matrix transpose() const
    {
        Matrix out;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                out(a,b) = data[b][a];
        return out;
    }

    __host__ __device__ double trace() const { return data[0][0] + data[1][1] + data[2][2]; }
    __host__ __device__ Matrix cwiseProduct(const Matrix &other) const
    {
        Matrix out;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                out(a,b) = data[a][b] * other(a,b);
        return out;
    }

    __host__ __device__ double sum() const
    {
        double result = 0;
        for (int b = 0; b < 3; ++b)
            for (int a = 0; a < 3; ++a)
                result += data[a][b];
        return result;
    }

    __host__ __device__ double squaredNorm() const { return cwiseProduct(*this).sum(); }
    __host__ __device__ double norm() const { return std::sqrt(squaredNorm()); }

    __host__ __device__ double cofactor(int a, int b) const
    {
        const int a1=(a+1)%3, a2=(a+2)%3, b1=(b+1)%3, b2=(b+2)%3;
        return data[a1][b1]*data[a2][b2] - data[a1][b2]*data[a2][b1];
    }

    __host__ __device__ double determinant() const
    {
        return data[0][0]*cofactor(0,0) + data[0][1]*cofactor(0,1) + data[0][2]*cofactor(0,2);
    }

    __host__ __device__ Matrix inverse() const
    {
        const double inverse_determinant = 1.0 / determinant();
        Matrix out;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                out(a,b) = cofactor(b,a)*inverse_determinant;
        return out;
    }
};
struct History { Matrix green, viscous; double dissipation; };
struct Response { Matrix pk1, cauchy; double free_energy; };

__device__ Matrix loadGelTensor(const ShapeHessian3 &input)
{
    Matrix value;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            value(a,b) = input.derivative[a][b];
    return value;
}

__device__ ShapeHessian3 storeGelTensor(const Matrix &value)
{
    ShapeHessian3 result;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            result.derivative[a][b] = value(a,b);
    return result;
}

__global__ void evaluateGelPoints(const GelPointState *points, std::size_t count, GelPointParameters parameters,
                                  double dt, GelPointResult *result)
{
    const auto i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= count)
        return;

    const auto &point = points[i];
    const Matrix F = loadGelTensor(point.deformation);
    const History old{loadGelTensor(point.green), loadGelTensor(point.viscous), point.dissipation};
    GelPointResult value{};
    History next;
    Response response;
    value.status = advanceGelKernel(F, dt, parameters, old, next);
    if (value.status == GelKernelStatus::OK)
        value.status = gelResponseKernel(F, parameters, next, response);
    if (value.status == GelKernelStatus::OK)
    {
        value.green = storeGelTensor(next.green);
        value.viscous = storeGelTensor(next.viscous);
        value.pk1 = storeGelTensor(response.pk1);
        value.cauchy = storeGelTensor(response.cauchy);
        value.dissipation = next.dissipation;
        value.free_energy = response.free_energy;
    }
    result[i] = value;
}

class CudaGelPoints final : public GelPointBatch
{
  public:
    CudaGelPoints()
    {
        int device = 0, limit = 0;
        check(cudaGetDevice(&device), "select gel material device");
        check(cudaDeviceGetAttribute(&limit, cudaDevAttrMaxGridDimX, device), "read gel material grid");
        max_blocks_ = static_cast<std::size_t>(limit);
    }

    void evaluate(const std::vector<GelPointState> &points, GelPointParameters parameters,
                  double dt, std::vector<GelPointResult> &result) override
    {
        const auto count = points.size();
        result.resize(count);
        if (count == 0)
            return;
        constexpr unsigned threads = 128;
        const auto blocks = (count - 1) / threads + 1;
        if (blocks > max_blocks_)
            throw std::length_error("Gel material batch exceeds the device grid");
        upload(points_, points, "gel material state");
        result_.reserve(count);
        evaluateGelPoints<<<static_cast<unsigned>(blocks), threads>>>(points_.get(), count, parameters, dt, result_.get());
        check(cudaGetLastError(), "evaluate gel material points");
        download(result, result_, "gel material response");
    }

  private:
    std::size_t max_blocks_ = 0;
    DeviceBuffer<GelPointState> points_;
    DeviceBuffer<GelPointResult> result_;
};
}

std::unique_ptr<GelPointBatch> makeCudaGelPoints()
{
    return std::make_unique<CudaGelPoints>();
}
}
