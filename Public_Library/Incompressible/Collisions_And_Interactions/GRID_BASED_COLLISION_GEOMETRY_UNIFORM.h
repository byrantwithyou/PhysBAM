//#####################################################################
// Copyright 2005-2009, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_BASED_COLLISION_GEOMETRY_UNIFORM
//#####################################################################
#ifndef __GRID_BASED_COLLISION_GEOMETRY_UNIFORM__
#define __GRID_BASED_COLLISION_GEOMETRY_UNIFORM__

#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_MAC_HELPER.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY.h>
namespace PhysBAM{

template<class TV> class COLLISION_GEOMETRY;
template<class TV> class RAY;

template <class TV>
class GRID_BASED_COLLISION_GEOMETRY_UNIFORM:public GRID_BASED_COLLISION_GEOMETRY<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<int,TV_INT> T_ARRAYS_INT;
    typedef ARRAY<int,FACE_INDEX<TV::m> > T_FACE_ARRAYS_INT;
    typedef LINEAR_INTERPOLATION_MAC_HELPER<TV> T_LINEAR_INTERPOLATION_MAC_HELPER;
public:

    typedef GRID_BASED_COLLISION_GEOMETRY<TV> BASE;
    using BASE::collision_geometry_collection;using BASE::grid;using BASE::collision_thickness;using BASE::objects_in_cell;using BASE::cell_neighbors_visible;
    using BASE::face_neighbors_visible;using BASE::Get_Body_Penetration;using BASE::occupied_blocks;using BASE::swept_occupied_blocks;using BASE::Is_Active;
    using BASE::Latest_Crossover;using BASE::Any_Simplex_Crossover;using BASE::Intersection_With_Any_Simplicial_Object;

    const ARRAY<bool,TV_INT>* outside_fluid;
    bool use_collision_face_neighbors;

    GRID_BASED_COLLISION_GEOMETRY_UNIFORM(GRID<TV>& grid_input);
    virtual ~GRID_BASED_COLLISION_GEOMETRY_UNIFORM();

    bool Occupied_Cell_Center(const TV_INT& cell_index) const
    {return occupied_blocks(cell_index);} // this will check the occupied block that is left,bottom,front of the current cell.

    bool Swept_Occupied_Cell_Center(const TV_INT& cell_index) const
    {return swept_occupied_blocks(cell_index);} // this will check the occupied block that is left,bottom,front of the current cell.

    bool Occupied_Face_Center(const FACE_INDEX<TV::m>& face) const
    {return occupied_blocks(grid.Face_Node_Index(face.axis,face.index,0));}

    bool Swept_Occupied_Face_Center(const FACE_INDEX<TV::m>& face) const
    {return swept_occupied_blocks(grid.Face_Node_Index(face.axis,face.index,0));}

    bool Latest_Cell_Crossover(const TV_INT& cell_index,const T dt) const
    {COLLISION_GEOMETRY_ID body_id;int aggregate_id;TV initial_hit_point,X=grid.Center(cell_index);
    return Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point);}

    bool Latest_Face_Crossover(const TV_INT& face_index,const int axis,const T dt) const
    {COLLISION_GEOMETRY_ID body_id;int aggregate_id;TV initial_hit_point,X=grid.Face(FACE_INDEX<TV::m>(axis,face_index));
    return Latest_Crossover(X,X,dt,body_id,aggregate_id,initial_hit_point);}

    bool Latest_Velocity_Crossover(const int axis,const TV_INT& face_index,const T dt,T& face_velocity) const;

    bool Latest_Velocity_Crossover(const int side,const int axis,const TV_INT& face_index,const T dt,T& face_velocity) const;

    bool Any_Crossover(const TV& start_X,const TV& end_X,const T dt) const
    {return Any_Simplex_Crossover(start_X,end_X,dt);} // TODO: use object ids

//#####################################################################
    void Initialize_Grids();
    virtual void Compute_Occupied_Blocks(const bool with_body_motion,const T extra_thickness,const T body_thickness_factor);
    void Compute_Grid_Visibility();
    void Compute_Psi_N(ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,FACE_INDEX<TV::m> >* face_velocities=0) const;
    void Compute_Psi_N_Zero_Velocity(ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,FACE_INDEX<TV::m> >* face_velocities=0) const;
    void Compute_Simplices_In_Cell(ARRAY<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >,TV_INT>& simplices_in_cell,const ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>& bodies,
        int ghost_cells,T thickness,bool assume_active) const;
    bool Inside_Any_Simplex_Of_Any_Body(const TV& location,COLLISION_GEOMETRY_ID& body_id,int& aggregate_id) const;
    bool Inside_Any_Body(const TV& location,COLLISION_GEOMETRY_ID& body_id) const;
    bool Implicit_Geometry_Lazy_Inside_Any_Body(const TV& location,COLLISION_GEOMETRY_ID& body_id) const;
    bool Closest_Non_Intersecting_Point_Of_Any_Body(RAY<TV>& ray,COLLISION_GEOMETRY_ID& body_id) const;
    bool Cell_Center_Visible_From_Face(const TV_INT& cell,const int axis,const TV_INT& face_index) const;
    bool Face_Velocity(const int axis,const TV_INT& face_index,const TV_INT* cells,const int number_of_cells,const TV& X,T& face_velocity) const;
    bool Face_Velocity(const int side,const int axis,const TV_INT& face_index,const TV_INT* cells,const int number_of_cells,const TV& X,T& face_velocity) const;
    bool Cell_Center_Intersection(const TV_INT& cell_index,const TV_INT* cell_indices_for_body_id,
        const int number_of_cells_for_body_id,const TV& X,COLLISION_GEOMETRY_ID& body_id,int& aggregate_id,
        TV& intersection_point) const;
    void Compute_Simplices_In_Cell(ARRAY<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >,TV_INT>& simplices_in_cell,
        int ghost_cells,T thickness) const;
    bool Latest_Cell_Crossover_And_Velocity(const TV_INT& cell_index,const T dt,TV& velocity) const;
    TV Object_Velocity(const COLLISION_GEOMETRY_ID body_id,const int aggregate_id,const TV& X) const;
//#####################################################################
};
}
#endif
