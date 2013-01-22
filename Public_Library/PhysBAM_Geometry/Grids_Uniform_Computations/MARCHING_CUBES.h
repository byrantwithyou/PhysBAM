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
template<class T> class POINT_SIMPLICES_1D;
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
    enum WORKAROUND {max_elements=(d==3?5:2),num_corners=1<<d,num_edges=d<<(d-1),num_pts=num_corners+num_edges};
    unsigned short elements[max_elements];

    static int edge_lookup[num_corners][num_corners];
    static int vertex_lookup[num_edges][2];
    static int permute_map[d==3?6:2][num_pts];

    static int face_map[2*d][(1<<(d-1))+((d-1)<<(d>1?d-2:0))];
    static const ARRAY<MARCHING_CUBES_CASE<d> >& Case_Table();
};

template<int d>
struct MARCHING_CUBES_INTERIOR_CASE
{
    enum WORKAROUND {max_elements=(d==3?16:(d==2?4:1)),num_corners=1<<d,num_edges=d<<(d-1),num_pts=num_corners+num_edges};
    unsigned int elements[2][max_elements];

    static const ARRAY<MARCHING_CUBES_INTERIOR_CASE<d> >& Case_Table();
};

template<class TV>
class MARCHING_CUBES
{
public:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
    enum WORKAROUND {num_corners=1<<TV::m,num_edges=TV::m<<(TV::m-1),num_pts=num_corners+num_edges};
    enum WORKAROUND1 {boundary_corners=1<<(TV::m-1),boundary_edges=(TV::m==3?4:1),boundary_pts=boundary_corners+boundary_edges};

    MARCHING_CUBES() {}
    ~MARCHING_CUBES() {}

    static void Compute_Phis_For_Cell(VECTOR<T,num_corners>& phis,const ARRAY<T,TV_INT>& phi,const TV_INT& cell);
    static int Compute_Points_For_Cell(VECTOR<TV,num_pts>& pts,const VECTOR<T,num_corners>& phis);
    static void Get_Elements_For_Cell(ARRAY<VECTOR<TV,TV::m> >& surface,const VECTOR<VECTOR<ARRAY<VECTOR<TV,TV::m> >*,2>,2*TV::m>& boundary,const VECTOR<T,num_corners>& phis);
    static void Get_Elements_For_Cell(ARRAY<VECTOR<TV,TV::m> >& surface,const VECTOR<VECTOR<ARRAY<VECTOR<TV,TV::m> >*,2>,2*TV::m>& boundary,const ARRAY<T,TV_INT>& phi,const TV_INT& cell);
    static void Get_Interior_Elements_For_Cell(ARRAY<VECTOR<TV,TV::m+1> >* interior,ARRAY<VECTOR<TV,TV::m+1> >* exterior,const VECTOR<T,num_corners>& phis);
    static void Get_Interior_Elements_For_Cell(ARRAY<VECTOR<TV,TV::m+1> >* interior,ARRAY<VECTOR<TV,TV::m+1> >* exterior,const ARRAY<T,TV_INT>& phi,const TV_INT& cell);
    static int Create_Surface(T_SURFACE& surface,const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi);
//#####################################################################
};
}
#endif
