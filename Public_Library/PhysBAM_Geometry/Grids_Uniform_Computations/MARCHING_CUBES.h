//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MARCHING_CUBES
//#####################################################################
#ifndef __MARCHING_CUBES__
#define __MARCHING_CUBES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class T> class TRIANGULATED_SURFACE;
template<class T> class SEGMENTED_CURVE;
template<class T> class TRIANGLE_3D;
template<class T> class SEGMENT_2D;
template<class TV> class GRID;

/*
  scccccbbbbbaaaaa
  aaaaa, bbbbb, ccccc = vertex index (12-19 cube vertices, 0-11 cube edges)
  s = start new sheet (0 or 1)
 */
template<int d>
struct MARCHING_CUBES_CASE
{
    enum WORKAROUND {max_surface=(d==3?5:2),max_boundary=(d==3?5:2),num_corners=1<<d,num_edges=d<<(d-1),num_pts=num_corners+num_edges};
    unsigned short surface[max_surface];
    unsigned short boundary[max_boundary];
    unsigned short proj_dir;
    unsigned short enclose_inside;

    static int edge_lookup[num_corners][num_corners];
    static int vertex_lookup[num_edges][2];
    static int permute_map[4*d-6][num_pts];
};

template<class TV>
class MARCHING_CUBES
{
public:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    enum WORKAROUND {num_corners=1<<TV::m,num_edges=TV::m<<(TV::m-1),num_pts=num_corners+num_edges};

    MARCHING_CUBES() {}
    ~MARCHING_CUBES() {}

    static const ARRAY<MARCHING_CUBES_CASE<TV::m> >& Case_Table();
    static void Initialize_Case_Table(ARRAY<MARCHING_CUBES_CASE<2> >& table);
    static void Initialize_Case_Table(ARRAY<MARCHING_CUBES_CASE<3> >& table);
    static void Initialize_Neighbor_Cases(ARRAY<MARCHING_CUBES_CASE<TV::m> >& table, int c);
    static void Get_Elements_For_Cell(ARRAY<T_FACE>& surface,ARRAY<T_FACE>& boundary,int& direction,bool& enclose_inside,
        const ARRAY<T,TV_INT>& phi,const TV_INT& cell);
    static int Create_Surface(T_SURFACE& surface,const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi);
//#####################################################################
};
}
#endif
