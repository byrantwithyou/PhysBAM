//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIS_STENCIL_BOUNDARY_UNIFORM
//#####################################################################
#ifndef __BASIS_STENCIL_BOUNDARY_UNIFORM__
#define __BASIS_STENCIL_BOUNDARY_UNIFORM__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Symbolics/MULTIVARIATE_POLYNOMIAL.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class TV>
struct BASIS_STENCIL_BOUNDARY_UNIFORM
{
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_OBJECT;

    // Assume constant stencil, variables indexed per element
    const T_OBJECT& object;

    BASIS_STENCIL_BOUNDARY_UNIFORM(const T_OBJECT& object_input): object(object_input) {}
    ~BASIS_STENCIL_BOUNDARY_UNIFORM() {}
};
}
#endif
