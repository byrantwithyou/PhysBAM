//#####################################################################
// Copyright 2006-2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLEX_OBJECT 
//#####################################################################
#ifndef __SIMPLEX_OBJECT__
#define __SIMPLEX_OBJECT__

#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{
template<class TV,int d> using
SIMPLEX_OBJECT=typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d>::OBJECT;
}
#endif
