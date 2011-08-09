//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PARSE_STRING_LEVEL_SET_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PARSE_STRING_LEVEL_SET_HPP

#include <string>

#include <boost/function.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
int Parse_String_Level_Set(
    const std::string& level_set_str,
    boost::function< T ( const VECTOR<T,D>& ) >& phi);

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_PARSING_PARSE_STRING_LEVEL_SET_HPP
