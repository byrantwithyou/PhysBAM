//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DENSITY_CONTAINER
//#####################################################################
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Incompressible/Advection_Collidable/ADVECTION_WRAPPER_COLLIDABLE_CELL.h>
#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <Incompressible/Grid_Based_Fields/DENSITY_CONTAINER.h>
using namespace PhysBAM;
template<class TV> DENSITY_CONTAINER<TV>::
DENSITY_CONTAINER(GRID<TV>& grid_input)
    :GRID_AND_ARRAY_CONTAINER<TV,T>(grid_input),density(array),nested_semi_lagrangian_collidable(0),semi_lagrangian_collidable(0)
{
    Set_Ambient_Density();
    boundary_default.Set_Fixed_Boundary(true,ambient_density);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DENSITY_CONTAINER<TV>::
~DENSITY_CONTAINER()
{
    delete nested_semi_lagrangian_collidable;delete semi_lagrangian_collidable;
}
//#####################################################################
// Function Use_Semi_Lagrangian_Collidable_Advection
//#####################################################################
template<class TV> void DENSITY_CONTAINER<TV>::
Use_Semi_Lagrangian_Collidable_Advection(const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list,const ARRAY<bool,FACE_INDEX<TV::m> >& face_velocities_valid_mask_input)
{
    assert(!nested_semi_lagrangian_collidable&&!semi_lagrangian_collidable);
    nested_semi_lagrangian_collidable=new T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL(body_list,valid_mask_current,valid_mask_next,ambient_density,false);
    semi_lagrangian_collidable=new ADVECTION_WRAPPER_COLLIDABLE_CELL<TV,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL,T_FACE_LOOKUP_COLLIDABLE>(*nested_semi_lagrangian_collidable,body_list,
        face_velocities_valid_mask_input);
    Set_Custom_Advection(*semi_lagrangian_collidable);
}
//#####################################################################
// Function Initialize_Array
//#####################################################################
template<class TV> void DENSITY_CONTAINER<TV>::
Initialize_Array(const int ghost_cells)
{
    GRID_AND_ARRAY_CONTAINER<TV,T>::Initialize_Array(ghost_cells);
    if(semi_lagrangian_collidable){
        valid_mask_current.Resize(grid.Cell_Indices(3),use_init,true);
        valid_mask_next.Resize(grid.Cell_Indices(3),no_init);}
}
//#####################################################################
// Function Fill_Beta_At_Faces
//#####################################################################
template<class TV> void DENSITY_CONTAINER<TV>::
Fill_Beta_At_Faces(const T dt,const T time,ARRAY<T,FACE_INDEX<TV::m> >& beta_face) const
{
    ARRAY<T,TV_INT> density_ghost(grid.Cell_Indices(1),no_init);
    boundary->Fill_Ghost_Cells(grid,density,density_ghost,dt,time,1);
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        TV_INT face_index=iterator.Face_Index();
        TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
        T density_face=(density_ghost(first_cell_index)+density_ghost(second_cell_index))*(T).5;
        beta_face.Component(axis)(face_index)=(T)1/density_face;}
}
//#####################################################################
// Function Get_Ghost_Density
//#####################################################################
template<class TV> void DENSITY_CONTAINER<TV>::
Get_Ghost_Density(const T dt,const T time,const int number_of_ghost_cells,ARRAY<T,TV_INT>& density_ghost) const
{
    boundary->Fill_Ghost_Cells(grid,density,density_ghost,dt,time,number_of_ghost_cells);
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class TV> void DENSITY_CONTAINER<TV>::
Euler_Step(const T dt,const T time,const int number_of_ghost_cells)
{  
    GRID_AND_ARRAY_CONTAINER<TV,T>::Euler_Step(dt,time,number_of_ghost_cells);
    array.Clamp_Below(0); // density needs to be non-negative
}
//#####################################################################
namespace PhysBAM{
template class DENSITY_CONTAINER<VECTOR<float,1> >;
template class DENSITY_CONTAINER<VECTOR<float,2> >;
template class DENSITY_CONTAINER<VECTOR<float,3> >;
template class DENSITY_CONTAINER<VECTOR<double,1> >;
template class DENSITY_CONTAINER<VECTOR<double,2> >;
template class DENSITY_CONTAINER<VECTOR<double,3> >;
}
