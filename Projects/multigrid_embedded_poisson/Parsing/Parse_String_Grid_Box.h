//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PARSE_STRING_GRID_BOX_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PARSE_STRING_GRID_BOX_HPP

#include <string>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
int Parse_String_Grid_Box(
    const std::string& grid_box_str,
    T (&min_x)[D], T (&max_x)[D]);

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PARSE_STRING_GRID_BOX_HPP
