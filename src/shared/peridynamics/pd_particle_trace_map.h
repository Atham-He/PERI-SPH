#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace SPH::peridynamics
{
struct ParticleTraceEntry
{
    std::size_t node;
    double weight;
};

struct ParticleTraceMap
{
    std::size_t node_count;
    std::vector<std::size_t> offsets;
    std::vector<ParticleTraceEntry> entries;
};

struct ParticleTraceContribution
{
    std::size_t row;
    double weight;
};

struct ParticleTraceTranspose
{
    std::vector<std::size_t> offsets;
    std::vector<ParticleTraceContribution> entries;
};

inline void validateParticleTrace(const ParticleTraceMap &map, std::size_t rows)
{
    if (map.offsets.size() != rows + 1 || map.offsets.front() != 0 ||
        map.offsets.back() != map.entries.size())
        throw std::invalid_argument("Invalid particle trace offsets");
    for (std::size_t i = 1; i < map.offsets.size(); ++i)
        if (map.offsets[i] < map.offsets[i - 1])
            throw std::invalid_argument("Unordered particle trace offsets");
    for (const auto &entry : map.entries)
        if (entry.node >= map.node_count)
            throw std::invalid_argument("Particle trace contains a node outside its port");
}

template <class ActiveRow>
ParticleTraceTranspose transposeParticleTrace(const ParticleTraceMap &map, ActiveRow active, double sign = 1.0)
{
    ParticleTraceTranspose result;
    result.offsets.assign(map.node_count + 1, 0);
    const std::size_t rows = map.offsets.size() - 1;
    for (std::size_t row = 0; row < rows; ++row)
        if (active(row))
            for (auto k = map.offsets[row]; k < map.offsets[row + 1]; ++k)
                ++result.offsets[map.entries[k].node + 1];
    for (std::size_t i = 1; i < result.offsets.size(); ++i)
        result.offsets[i] += result.offsets[i - 1];

    result.entries.resize(result.offsets.back());
    auto cursor = result.offsets;
    for (std::size_t row = 0; row < rows; ++row)
        if (active(row))
            for (auto k = map.offsets[row]; k < map.offsets[row + 1]; ++k)
            {
                const auto entry = map.entries[k];
                result.entries[cursor[entry.node]++] = {row, sign * entry.weight};
            }
    return result;
}
}
