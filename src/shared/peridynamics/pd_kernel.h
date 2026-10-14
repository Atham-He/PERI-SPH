/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	pd_kernel.h
 * @brief 	Influence function (kernel) for peridynamics.
 * @details	Ported from Icebreaking_PD-SPH (Haotian Shi, Luwen Zhang, SJTU)
 * 			to the present SPHinXsys variable-registration framework.
 */

#ifndef PD_KERNEL_H
#define PD_KERNEL_H

#include "base_kernel.h"

namespace SPH
{
/**
 * @class Winfunc
 * @brief Influence function for peridynamics.
 * @details Cubic B-spline shaped weight with support 2h (the horizon),
 * normalization follows the monograph of Silling & Askari (2005).
 */
class Winfunc : public Kernel
{
  public:
    explicit Winfunc(Real h);
    virtual ~Winfunc() {};

    virtual Real W_1D(const Real q) const override;
    virtual Real W_2D(const Real q) const override;
    virtual Real W_3D(const Real q) const override;

    virtual Real dW_1D(const Real q) const override;
    virtual Real dW_2D(const Real q) const override;
    virtual Real dW_3D(const Real q) const override;

    virtual Real d2W_1D(const Real q) const override;
    virtual Real d2W_2D(const Real q) const override;
    virtual Real d2W_3D(const Real q) const override;
};
} // namespace SPH
#endif // PD_KERNEL_H
