//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARAMS_GENERAL_PARAMS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARAMS_GENERAL_PARAMS_HPP

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

struct GENERAL_PARAMS
{
    unsigned int rng_seed;
    unsigned int n_thread;
    bool randomized_check;

    GENERAL_PARAMS()
        : rng_seed(5489),
          n_thread(1),
          randomized_check(false)
    { }
};

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARAMS_GENERAL_PARAMS_HPP
