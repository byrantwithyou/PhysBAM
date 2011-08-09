//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARAMS_GRID_PARAMS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARAMS_GRID_PARAMS_HPP

#include <algorithm>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
struct GRID_PARAMS
{
    T min_x[D];
    T max_x[D];
    unsigned int n_cell[D];

    GRID_PARAMS()
    {
        std::fill(&min_x[0], &min_x[D], static_cast<T>(-1));
        std::fill(&max_x[0], &max_x[D], static_cast<T>(+1));
        std::fill(&n_cell[0], &n_cell[D], static_cast< unsigned int >(1));
    }
};

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARAMS_GRID_PARAMS_HPP
