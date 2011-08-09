//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PARSE_STRING_SOLVER_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PARSE_STRING_SOLVER_HPP

#include <string>

#include "../Params/SOLVER_PARAMS.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

int Parse_String_Solver(
    const std::string& solver_str,
    SOLVER_PARAMS& solver);

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PARSE_STRING_SOLVER_HPP
