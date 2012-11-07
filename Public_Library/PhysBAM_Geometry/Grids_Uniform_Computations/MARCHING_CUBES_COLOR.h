//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MARCHING_CUBES_COLOR
//#####################################################################
#ifndef __MARCHING_CUBES_COLOR__
#define __MARCHING_CUBES_COLOR__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
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

    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
    enum WORKAROUND {num_corners=1<<TV::m,num_edges=TV::m<<(TV::m-1),num_pts=num_corners+num_edges};

    MARCHING_CUBES_COLOR() {}
    ~MARCHING_CUBES_COLOR() {}

    static void Initialize_Case_Table();
    static void Get_Elements_For_Cell(ARRAY<TRIPLE<T_FACE,int,int> >& interface,ARRAY<PAIR<T_FACE,int> >& boundary,
        const VECTOR<int,num_corners>& colors,const VECTOR<T,num_corners>& phi);

//#####################################################################

    typedef HASHTABLE<int,T_SURFACE*> HASH_BOUNDARY;
    typedef HASHTABLE<int,INTERVAL<int> > HASH_CELL_BOUNDARY;
    typedef HASHTABLE<VECTOR<int,2>,T_SURFACE*> HASH_INTERFACE;
    typedef HASHTABLE<VECTOR<int,2>,INTERVAL<int> > HASH_CELL_INTERFACE;
    typedef HASHTABLE<TV_INT,PAIR<HASH_CELL_INTERFACE,HASH_CELL_BOUNDARY> > HASH_CELL_TO_ELEMENT;

    static void Get_Elements(const GRID<TV>& grid,HASH_INTERFACE& interface,HASH_BOUNDARY& boundary,
        HASH_CELL_TO_ELEMENT& cell_to_element,const ARRAY<int,TV_INT>& phi_color,const ARRAY<T,TV_INT>& phi_value,
        const int newton_steps=20,const bool verbose=false);

//#####################################################################
};
}
#endif
;
