//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARSING_PARSE_STRING_EXAMPLE_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARSING_PARSE_STRING_EXAMPLE_HPP

#include <string>

#include "../Params/EXAMPLE_PARAMS.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template< class T, int D >
int Parse_String_Example(
    const std::string& example_str,
    typename EXAMPLE_PARAMS<T,D>::PROBLEM_TYPE& problem);

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARSING_PARSE_STRING_EXAMPLE_HPP
