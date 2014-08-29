//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MARCHING_CUBES_COLOR
//#####################################################################
#ifndef __MARCHING_CUBES_COLOR__
#define __MARCHING_CUBES_COLOR__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>

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

    MARCHING_CUBES_COLOR(){}
    ~MARCHING_CUBES_COLOR(){}

    struct INTERFACE_ELEMENT
    {
        T_FACE face;
        VECTOR<int,2> color_pair;
    };

    struct BOUNDARY_ELEMENT
    {
        T_FACE face;
        int color;
    };

    struct CELL_ELEMENTS
    {
        ARRAY<INTERFACE_ELEMENT> interface;
        ARRAY<BOUNDARY_ELEMENT> boundary;
    };

    // Boundary elements have standard orientation
    //   2D: regions are traced counterclockwise
    //   3D: triangles are counterclockwise when viewed from the outside
    // Interface elements are consistent with boundary elements from first color
    static void Initialize_Case_Table();
    static void Get_Elements_For_Cell(ARRAY<INTERFACE_ELEMENT>& interface,ARRAY<BOUNDARY_ELEMENT>& boundary,
        const VECTOR<int,num_corners>& colors,const VECTOR<T,num_corners>& phi);
    static void Get_Elements_For_Cell(ARRAY<INTERFACE_ELEMENT>& interface,ARRAY<BOUNDARY_ELEMENT>& boundary,
        const VECTOR<int,num_corners>& colors,const VECTOR<T,num_corners>& phi,const RANGE<TV>& cell_range);

//#####################################################################

    typedef HASHTABLE<int,T_SURFACE*> HASH_BOUNDARY;
    typedef HASHTABLE<int,INTERVAL<int> > HASH_CELL_BOUNDARY;
    typedef HASHTABLE<VECTOR<int,2>,T_SURFACE*> HASH_INTERFACE;
    typedef HASHTABLE<VECTOR<int,2>,INTERVAL<int> > HASH_CELL_INTERFACE;

    struct HASH_CELL_DATA
    {
        HASH_CELL_INTERFACE interface;
        VECTOR<HASH_CELL_BOUNDARY,TV::m> boundary;
    };

//#####################################################################
    static void Get_Elements(HASHTABLE<TV_INT,CELL_ELEMENTS>& index_to_cell_elements,const GRID<TV>& grid,
        const ARRAY<int,TV_INT>& phi_color,const ARRAY<T,TV_INT>& phi_value,
        const int newton_steps=20,const bool verbose=false);

//#####################################################################

private:

    static void Fix_Mesh(GEOMETRY_PARTICLES<TV>& particles,ARRAY<int>& particle_dofs,HASHTABLE<TV_INT>& variable_cells,
        const HASH_INTERFACE& interface,const HASH_BOUNDARY& boundary,const HASHTABLE<TV_INT,HASH_CELL_DATA>& index_to_cell_data,
        const HASHTABLE<FACE_INDEX<TV::m>,int>& edge_vertices,const HASHTABLE<FACE_INDEX<TV::m>,int>& face_vertices,
        const HASHTABLE<TV_INT,int>& cell_vertices,const HASHTABLE<TV_INT,int>& node_vertices,
        const HASHTABLE<TV_INT>& junction_cells,const int iterations,const bool verbose);

    static void Save_Mesh(HASHTABLE<TV_INT,CELL_ELEMENTS>& index_to_cell_elements,const GRID<TV>& grid,const HASH_INTERFACE& interface,
        const HASH_BOUNDARY& boundary,const HASHTABLE<TV_INT,HASH_CELL_DATA>& index_to_cell_data,const GEOMETRY_PARTICLES<TV>& particles,
        const bool recut_cells=false,const ARRAY<int>* const particle_dofs=0,const HASHTABLE<TV_INT>* const variable_cells=0);

//#####################################################################
};
}
#endif
