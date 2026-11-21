#pragma once

#include "pd_crack_strength_batch.h"
#include <cmath>

#if defined(__CUDACC__)
#define SPH_CRACK_STRENGTH_HD __host__ __device__
#else
#define SPH_CRACK_STRENGTH_HD
#endif

namespace SPH::peridynamics
{
/** Strength screening of a sampled material plane in the current configuration. */
SPH_CRACK_STRENGTH_HD inline CrackStrengthResult crackStrength(
    const CrackStrengthState &first, const CrackStrengthState &second,
    const double (&reference_normal)[3], double strength)
{
    double deformation[3][3], traction[3]{};
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
        {
            deformation[a][b] = 0.5 * (first.deformation[a][b] + second.deformation[a][b]);
            traction[a] += 0.5 * (first.piola[a][b] + second.piola[a][b]) * reference_normal[b];
        }

    // Nanson's area vector, cof(F) * n0. Check the determinant explicitly,
    // since a singular deformation has no admissible current material plane.
    double normal[3]{};
    double determinant = 0;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
        {
            const int row1 = (a + 1) % 3, row2 = (a + 2) % 3;
            const int col1 = (b + 1) % 3, col2 = (b + 2) % 3;
            const double cofactor = deformation[row1][col1] * deformation[row2][col2] -
                                    deformation[row1][col2] * deformation[row2][col1];
            normal[a] += cofactor * reference_normal[b];
            if (a == 0)
                determinant += deformation[0][b] * cofactor;
        }
    if (!std::isfinite(determinant) || determinant == 0 ||
        !std::isfinite(normal[0]) || !std::isfinite(normal[1]) || !std::isfinite(normal[2]))
        return {0, false};

    const double norm = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    if (norm == 0)
        return {0, false};

    double opening_traction = 0, traction_squared = 0;
    for (int a = 0; a < 3; ++a)
    {
        opening_traction += traction[a] * (normal[a] / norm);
        traction_squared += traction[a] * traction[a];
    }
    return {opening_traction > 0 ? std::sqrt(traction_squared) / strength : 0, true};
}
}

#undef SPH_CRACK_STRENGTH_HD
