//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PARSE_STRING_N_CELL_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PARSE_STRING_N_CELL_HPP

#include <string>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< int D >
int Parse_String_N_Cell(
    const std::string& n_cell_str,
    unsigned int (&n_cell)[D]);

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PARSE_STRING_N_CELL_HPP
