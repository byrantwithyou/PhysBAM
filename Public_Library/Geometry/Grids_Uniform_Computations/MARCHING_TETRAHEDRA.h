//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MARCHING_TETRAHEDRA
//#####################################################################
#ifndef __MARCHING_TETRAHEDRA__
#define __MARCHING_TETRAHEDRA__

#include <Core/Arrays/ARRAY.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class T> class TRIANGULATED_SURFACE;
template<class T> class SEGMENTED_CURVE;
template<class T> class TRIANGLE_3D;
template<class T> class SEGMENT_2D;
template<class TV> class GRID;

// 3d: 0-3 = vectex, 4=01 5=02 6=03 7=12 8=13 9=23
// 2d: 0-2 = vectex, 3=01 4=02 5=12
// ccccbbbbaaaa
template<int d>
struct MARCHING_TETRAHEDRA_CASE
{
    enum WORKAROUND {max_surface=(d==3?2:1),max_boundary=(d==3?2:1),num_corners=d+1,num_edges=d*(d+1)/2,num_pts=num_corners+num_edges};
    unsigned short surface[max_surface];
    unsigned short boundary[d+1][2][max_boundary]; // [face][in/out][elements]

    static int vertex_lookup[num_edges][2];
};

template<class TV>
class MARCHING_TETRAHEDRA
{
public:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    enum WORKAROUND {num_corners=TV::m+1,num_edges=TV::m*(TV::m+1)/2,num_pts=num_corners+num_edges};

    MARCHING_TETRAHEDRA() {}
    ~MARCHING_TETRAHEDRA() {}

    static const MARCHING_TETRAHEDRA_CASE<TV::m>* Case_Table();
    static void Get_Elements_For_Tetrahedron(ARRAY<T_FACE>& surface,const VECTOR<VECTOR<ARRAY<T_FACE>*,2>,TV::m+1>& boundary,const VECTOR<T,TV::m+1>& phi,const VECTOR<TV,TV::m+1>& X);
    static void Fill_Faces(ARRAY<T_FACE>& faces,const unsigned short* face_encoding,const TV* pts);
//#####################################################################
};
}
#endif
