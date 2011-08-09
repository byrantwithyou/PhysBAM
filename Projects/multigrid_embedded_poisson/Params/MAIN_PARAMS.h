//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARAMS_MAIN_PARAMS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARAMS_MAIN_PARAMS_HPP

#include "EXAMPLE_PARAMS.h"
#include "GENERAL_PARAMS.h"
#include "GRID_PARAMS.h"
#include "LEVEL_SET_PARAMS.h"
#include "OUTPUT_PARAMS.h"
#include "SOLVER_PARAMS.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
struct MAIN_PARAMS
{
    GENERAL_PARAMS general;
    GRID_PARAMS<T,D> grid;
    LEVEL_SET_PARAMS<T,D> level_set;
    EXAMPLE_PARAMS<T,D> example;
    SOLVER_PARAMS solver;
    OUTPUT_PARAMS output;
};

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARAMS_MAIN_PARAMS_HPP
