//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_MAKER_UNIFORM
//##################################################################### 
#ifndef __LEVELSET_MAKER_UNIFORM__
#define __LEVELSET_MAKER_UNIFORM__

#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class TV> class GRID;
template<class T> class SEGMENTED_CURVE_2D;
template<class T> class TRIANGULATED_SURFACE;

template<class TV>
class LEVELSET_MAKER_UNIFORM
{
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
public:
//#####################################################################
    static void Compute_Level_Set(T_SURFACE& surface,GRID<TV>& grid,int ghost_cells,ARRAY<T,TV_INT>& phi);
//#####################################################################
};
}
#endif
