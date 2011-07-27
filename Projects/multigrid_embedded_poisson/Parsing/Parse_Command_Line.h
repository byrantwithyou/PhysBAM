//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARSING_PARSE_COMMAND_LINE_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARSING_PARSE_COMMAND_LINE_HPP

#include "../Params/MAIN_PARAMS.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template< class T, int D >
int Parse_Command_Line(
    int argc, char* argv[],
    MAIN_PARAMS<T,D>& main_params);

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARSING_PARSE_COMMAND_LINE_HPP
