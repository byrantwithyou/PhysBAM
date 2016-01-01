//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SIGNED_DISTANCE
//##################################################################### 
#ifndef __TRIANGULATED_SURFACE_SIGNED_DISTANCE_UNIFORM__
#define __TRIANGULATED_SURFACE_SIGNED_DISTANCE_UNIFORM__
#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
namespace PhysBAM{
template<class TV> class GRID;

namespace SIGNED_DISTANCE{
template<class T> void Calculate(TRIANGULATED_SURFACE<T>& surface,const GRID<VECTOR<T,3> >& grid,ARRAY<T,VECTOR<int,3> >& phi,bool print_progress=false);
};
};
#endif
