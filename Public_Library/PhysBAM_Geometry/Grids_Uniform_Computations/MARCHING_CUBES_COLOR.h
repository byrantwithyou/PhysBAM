//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MARCHING_CUBES_COLOR
//#####################################################################
#ifndef __MARCHING_CUBES_COLOR__
#define __MARCHING_CUBES_COLOR__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class T> class TRIANGLE_3D;
template<class T> class SEGMENT_2D;
template<class TV> class GRID;

template<class TV>
class MARCHING_CUBES_COLOR
{
public:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
    enum WORKAROUND {num_corners=1<<TV::m,num_edges=TV::m<<(TV::m-1),num_pts=num_corners+num_edges};

    MARCHING_CUBES_COLOR() {}
    ~MARCHING_CUBES_COLOR() {}

    static void Initialize_Case_Table();
    static void Get_Elements_For_Cell(ARRAY<TRIPLE<T_FACE,int,int> >& surface,ARRAY<PAIR<T_FACE,int> >& boundary,
        const VECTOR<int,num_corners>& colors,const VECTOR<T,num_corners>& phi);
    static void Get_Elements(const GRID<TV>& grid,HASHTABLE<VECTOR<int,2>,T_SURFACE*>& surface,HASHTABLE<int,T_SURFACE*>& boundary,
        const ARRAY<int,TV_INT>& color,const ARRAY<T,TV_INT>& phi);
//#####################################################################
};
}
#endif
