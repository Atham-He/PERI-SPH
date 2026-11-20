#pragma once

#include <cuda_runtime.h>
#include <stdexcept>
#include <string>

namespace SPH::peridynamics::cuda_detail
{
inline void check(cudaError_t status, const char *operation)
{
    if (status != cudaSuccess)
        throw std::runtime_error(std::string("PERI-SPH CUDA ") + operation +
                                 " failed: " + cudaGetErrorString(status));
}

template <class T>
class DeviceBuffer
{
  public:
    ~DeviceBuffer() { cudaFree(data_); }
    DeviceBuffer() = default;
    DeviceBuffer(const DeviceBuffer &) = delete;
    DeviceBuffer &operator=(const DeviceBuffer &) = delete;

    T *get() const { return data_; }
    void reserve(std::size_t count)
    {
        if (count <= capacity_)
            return;
        check(cudaFree(data_), "buffer release");
        data_ = nullptr;
        capacity_ = 0;
        check(cudaMalloc(&data_, count * sizeof(T)), "buffer allocation");
        capacity_ = count;
    }

  private:
    T *data_ = nullptr;
    std::size_t capacity_ = 0;
};

}
