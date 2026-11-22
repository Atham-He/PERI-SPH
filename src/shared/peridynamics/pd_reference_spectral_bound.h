#pragma once

#include "data_type.h"
#include <Eigen/SparseCore>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace SPH::peridynamics
{
struct ReferenceSpectralBound
{
    Real upper = 0;
    Real largest_diagonal = 0;
};

/** Bound the mass-normalized gradient Gram G^T V G. The same scalar
 * spectrum bounds the three displacement components of the vector gradient.
 * Each input row gives neighbor coefficients; its center coefficient is their
 * negative sum, as in the pointwise differences used by the solid operator. */
template <class Rows>
ReferenceSpectralBound boundReferenceGradient(const Rows &rows, const Real *volume, const Real *mass)
{
    using Wide = long double;
    using Sparse = Eigen::SparseMatrix<Wide, Eigen::RowMajor>;
    using Index = Sparse::StorageIndex;
    using Triplet = Eigen::Triplet<Wide, Index>;
    const auto count = rows.size();
    if (count == 0 || count > std::size_t(std::numeric_limits<Index>::max()) / Dimensions)
        throw std::invalid_argument("Reference gradient spectral bound has invalid dimensions");

    std::vector<Wide> inverse_mass_root(count), error_rows(count, 0);
    std::vector<Triplet> triplets;
    std::size_t stencil_entries = count;
    for (std::size_t i = 0; i < count; ++i)
    {
        if (!std::isfinite(volume[i]) || volume[i] <= 0 || !std::isfinite(mass[i]) || mass[i] <= 0)
            throw std::invalid_argument("Reference gradient spectral bound requires positive volume and mass");
        inverse_mass_root[i] = 1 / std::sqrt(Wide(mass[i]));
        stencil_entries += rows[i].size();
    }
    triplets.reserve(Dimensions * stencil_entries);

    for (std::size_t i = 0; i < count; ++i)
    {
        Wide self[Dimensions]{};
        const Wide root_volume = std::sqrt(Wide(volume[i]));
        Wide self_enclosure = 0, neighbor_enclosure = 0;
        std::vector<Wide> neighbor_sizes;
        neighbor_sizes.reserve(rows[i].size());
        for (const auto &entry : rows[i])
        {
            if (entry.j >= count || entry.j == i || !entry.gradient.allFinite())
                throw std::invalid_argument("Invalid reference gradient neighbor");
            Wide size = 0;
            for (int axis = 0; axis < Dimensions; ++axis)
            {
                const Wide coefficient = entry.gradient[axis];
                self[axis] -= coefficient;
                size += std::abs(coefficient);
                if (coefficient != 0)
                    triplets.emplace_back(Index(Dimensions * i + axis), Index(entry.j),
                                          root_volume * coefficient * inverse_mass_root[entry.j]);
            }
            self_enclosure += size;
            neighbor_sizes.push_back(root_volume * size * inverse_mass_root[entry.j]);
            neighbor_enclosure += neighbor_sizes.back();
        }
        self_enclosure *= root_volume * inverse_mass_root[i];
        const Wide enclosure = self_enclosure + neighbor_enclosure;
        error_rows[i] += self_enclosure * enclosure;
        for (std::size_t j = 0; j < rows[i].size(); ++j)
            error_rows[rows[i][j].j] += neighbor_sizes[j] * enclosure;
        for (int axis = 0; axis < Dimensions; ++axis)
            if (self[axis] != 0)
                triplets.emplace_back(Index(Dimensions * i + axis), Index(i),
                                      root_volume * self[axis] * inverse_mass_root[i]);
    }

    Sparse gradient(Index(Dimensions * count), Index(count));
    gradient.setFromTriplets(triplets.begin(), triplets.end());
    Sparse gram = gradient.transpose() * gradient;

    // A similarity scaling cannot change eigenvalues. Each positive scale
    // yields a separate Gershgorin bound from the signed assembled Gram.
    // Iteration only improves the chosen scaling, with no convergence claim.
    constexpr Wide minimum_scale = 1.L / 32;
    std::vector<Wide> scale(count, 1), product(count);
    Wide largest_diagonal = 0;
    Wide upper = std::numeric_limits<Wide>::infinity();
    for (unsigned iteration = 0; iteration <= 12; ++iteration)
    {
        Wide largest = 0, scaled_upper = 0;
        for (Index i = 0; i < Index(count); ++i)
        {
            Wide sum = 0;
            for (typename Sparse::InnerIterator entry(gram, i); entry; ++entry)
            {
                if (iteration == 0 && entry.col() == i)
                    largest_diagonal = std::max(largest_diagonal, entry.value());
                sum += std::abs(entry.value()) * scale[entry.col()];
            }
            product[i] = sum;
            largest = std::max(largest, sum);
            scaled_upper = std::max(scaled_upper, sum / scale[i]);
        }
        if (!(largest > 0) || !std::isfinite(largest))
            throw std::runtime_error("Invalid reference gradient Gram bound");
        upper = std::min(upper, scaled_upper);
        for (std::size_t i = 0; i < count; ++i)
            scale[i] = std::max(minimum_scale, product[i] / largest);
    }

    // Wide accumulation plus an absolute rounding enclosure for coefficient
    // sums, products, and the scaled row sums. The l1 envelope also covers
    // cancellation in the implicit center coefficient.
    const Wide row_enclosure = *std::max_element(error_rows.begin(), error_rows.end());
    const Wide roundoff = 128 * std::numeric_limits<Wide>::epsilon() * (Dimensions * count + 1);
    if (!std::isfinite(row_enclosure) || roundoff >= .01)
        throw std::runtime_error("Reference gradient roundoff enclosure is unresolved");
    upper += roundoff * row_enclosure / minimum_scale;
    const Real rounded_upper = std::nextafter(Real(upper), std::numeric_limits<Real>::infinity());
    if (!(rounded_upper > 0) || !std::isfinite(rounded_upper))
        throw std::runtime_error("Reference gradient spectral bound is nonfinite");
    return {rounded_upper, Real(largest_diagonal)};
}
}
