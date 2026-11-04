/* Resource policy for the small gel control. No numerical threshold is relaxed. */
#ifndef PERI_GEL_EXECUTION_LIMITS_H
#define PERI_GEL_EXECUTION_LIMITS_H

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace gel_compression
{
class ResourceStop : public std::runtime_error
{
  public:
    ResourceStop(const std::string &reason, int code = 125)
        : std::runtime_error(reason), return_code(code) {}
    int return_code;
};

/** Prospective budget for the new <=8192-particle solid gel spatial study only.
 * 912 particles / 920 steps measured 2.764s and a 9.94MiB checkpoint; at fixed
 * horizon, doubling resolution predicts ~128x work and ~64x neighbor storage.
 * This does not modify the sealed R8 N<=1000 guard or any accuracy threshold.
 * Host storage must be identified separately from free space inside a WSL VHD. */
class ExecutionLimits
{
  public:
    static constexpr size_t max_particles = 8192;
    static constexpr size_t max_steps = 98304;
    static constexpr size_t max_threads = 8;

    ExecutionLimits() : start_(std::chrono::steady_clock::now()), checked_(start_)
    {
        const char *specified = std::getenv("PERI_HOST_STORAGE");
        host_storage_ = specified ? specified : "/mnt/e";
        checkStorage();
    }

    void check(size_t steps, double previous_wall_seconds)
    {
        if (steps >= max_steps)
            throw ResourceStop("accepted_step_budget");
        const auto now = std::chrono::steady_clock::now();
        if (previous_wall_seconds + elapsed() >= 3600)
            throw ResourceStop("wall_time_budget", 124);
        if (std::chrono::duration<double>(now - checked_).count() >= 5)
        {
            checkStorage();
            std::ifstream status("/proc/self/status");
            if (!status)
                throw ResourceStop("cannot_verify_process_memory");
            std::string line;
            bool found = false;
            while (std::getline(status, line))
                if (line.rfind("VmRSS:", 0) == 0)
                {
                    found = true;
                    if (std::stoull(line.substr(6)) >= 8ULL * 1024 * 1024)
                        throw ResourceStop("process_memory_budget");
                }
            if (!found)
                throw ResourceStop("cannot_verify_process_memory");
            checked_ = now;
        }
    }

    double elapsed() const
    {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
    }

  private:
    void checkStorage() const
    {
        constexpr uintmax_t gib = 1024ULL * 1024 * 1024;
        if (!std::filesystem::exists(host_storage_))
            throw ResourceStop("host_storage_unknown: set PERI_HOST_STORAGE to the actual host volume");
        if (std::filesystem::space(host_storage_).available < 120 * gib)
            throw ResourceStop("host_disk_reserve");
        if (std::filesystem::space(std::filesystem::current_path()).available < 150 * gib)
            throw ResourceStop("workspace_disk_reserve");
    }

    std::filesystem::path host_storage_;
    std::chrono::steady_clock::time_point start_, checked_;
};
} // namespace gel_compression
#endif
