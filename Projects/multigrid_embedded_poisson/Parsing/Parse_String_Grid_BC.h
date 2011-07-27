//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARSING_PARSE_STRING_GRID_BC_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARSING_PARSE_STRING_GRID_BC_HPP

#include <string>

#include "../Params/EXAMPLE_PARAMS.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

int Parse_String_Grid_BC(
    const std::string& grid_bc_str,
    EXAMPLE_PARAMS_BASE::GRID_BC_ID& grid_bc_id);

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARSING_PARSE_STRING_GRID_BC_HPP
