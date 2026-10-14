/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * Licensed under the Apache License, Version 2.0 (the "License");           *
 * you may not use this file except in compliance with the License.          *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	pd_body.h
 * @brief 	Peridynamic body and its adaptation.
 * @details	The PD body uses the Winfunc influence function and a horizon
 * 			of 1.5075 * spacing (support 2h = 3.015 * spacing), reproducing
 * 			the setup of Icebreaking_PD-SPH (Haotian Shi, Luwen Zhang, SJTU).
 */

#ifndef PD_BODY_H
#define PD_BODY_H

#include "adaptation.h"
#include "base_body.h"
#include "pd_kernel.h"

namespace SPH
{
/**
 * @class PDAdaptation
 * @brief Single-resolution adaptation with the peridynamic influence function.
 */
class PDAdaptation : public SPHAdaptation
{
  public:
    explicit PDAdaptation(Real global_resolution, Real h_spacing_ratio = 1.5075,
                          Real refinement_to_global = 1.0)
        : SPHAdaptation(global_resolution, h_spacing_ratio, refinement_to_global)
    {
        resetKernel<Winfunc>();
    };
    virtual ~PDAdaptation() {};
};

/**
 * @class PDBody
 * @brief A real body discretized by peridynamic particles.
 */
class PDBody : public RealBody
{
  public:
    template <typename... Args>
    explicit PDBody(Args &&...args)
        : RealBody(std::forward<Args>(args)...)
    {
        defineAdaptation<PDAdaptation>(1.5075);
    };
    virtual ~PDBody() {};
};
} // namespace SPH
#endif // PD_BODY_H
