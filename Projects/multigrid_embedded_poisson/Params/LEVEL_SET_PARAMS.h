//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARAMS_LEVEL_SET_PARAMS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARAMS_LEVEL_SET_PARAMS_HPP

#include <boost/function.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template< class T, int D >
struct LEVEL_SET_PARAMS
{
    float min_dist_to_vertex;
    boost::function< T ( const VECTOR<T,D>& ) > phi;

    LEVEL_SET_PARAMS()
        : min_dist_to_vertex(0)
    { }
};

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_PARAMS_LEVEL_SET_PARAMS_HPP
