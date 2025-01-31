//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TEMPERATURE_CONTAINER
//#####################################################################
#include <Core/Math_Tools/cube.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Incompressible/Advection_Collidable/ADVECTION_WRAPPER_COLLIDABLE_CELL.h>
#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <Incompressible/Grid_Based_Fields/TEMPERATURE_CONTAINER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TEMPERATURE_CONTAINER<TV>::
TEMPERATURE_CONTAINER(GRID<TV>& grid_input)
    :GRID_AND_ARRAY_CONTAINER<TV,T>(grid_input),temperature(array),nested_semi_lagrangian_collidable(0),semi_lagrangian_collidable(0)
{
    Set_Ambient_Temperature();
    boundary_default.Set_Fixed_Boundary(true,ambient_temperature);
    Set_Cooling_Constant();
    Set_Hot_Point();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> TEMPERATURE_CONTAINER<TV>::
~TEMPERATURE_CONTAINER()
{
    delete nested_semi_lagrangian_collidable;delete semi_lagrangian_collidable;
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class TV> void TEMPERATURE_CONTAINER<TV>::
Euler_Step(const T dt,const T time,const int number_of_ghost_cells)
{
    GRID_AND_ARRAY_CONTAINER<TV,T>::Euler_Step(dt,time,number_of_ghost_cells);
    array.Clamp_Below(ambient_temperature); // temperature needs to be above ambient temperature
    Apply_Cooling(dt,time);
}
//#####################################################################
// Function Apply_Cooling
//#####################################################################
template<class TV> void TEMPERATURE_CONTAINER<TV>::
Apply_Cooling(const T dt,const T time)
{  
    if(!cooling_constant) return;
    T constant=3*cooling_constant*dt/sqr(sqr(hot_point-ambient_temperature));
    for(CELL_ITERATOR<TV> iterator(grid,0);iterator.Valid();iterator.Next()) Apply_Individual_Cooling(temperature(iterator.Cell_Index()),constant);
}
// TODO: the following will go away once we figure out a way to merge all the iterators
template<class T,class TV> static void Apply_Cooling_Helper(TEMPERATURE_CONTAINER<TV>& container,const T dt,const T time)
{
    if(!container.cooling_constant) return;
    T constant=3*container.cooling_constant*dt/sqr(sqr(container.hot_point-container.ambient_temperature));
    for(int i=0;i<container.grid.number_of_cells;i++)container.Apply_Individual_Cooling(container.temperature(i),constant);
}
//#####################################################################
// Function Apply_Individual_Cooling
//#####################################################################
template<class TV> void TEMPERATURE_CONTAINER<TV>::
Apply_Individual_Cooling(T& temperature,const T constant)
{
    if(abs(temperature-ambient_temperature) > (T).1) temperature=ambient_temperature+pow(1/(constant+1/cube(temperature-ambient_temperature)),((T)1/3));
}
//#####################################################################
// Function Use_Semi_Lagrangian_Collidable_Advection
//#####################################################################
template<class TV> void TEMPERATURE_CONTAINER<TV>::
Use_Semi_Lagrangian_Collidable_Advection(const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list,const ARRAY<bool,FACE_INDEX<TV::m> >& face_velocities_valid_mask_input)
{
    assert(!nested_semi_lagrangian_collidable&&!semi_lagrangian_collidable);
    nested_semi_lagrangian_collidable=new T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL(body_list,valid_mask_current,valid_mask_next,ambient_temperature,false);
    semi_lagrangian_collidable=new ADVECTION_WRAPPER_COLLIDABLE_CELL<TV,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL,T_FACE_LOOKUP_COLLIDABLE>(*nested_semi_lagrangian_collidable,body_list,
        face_velocities_valid_mask_input);
    Set_Custom_Advection(*semi_lagrangian_collidable);
}
//#####################################################################
template class TEMPERATURE_CONTAINER<VECTOR<float,1> >;
template class TEMPERATURE_CONTAINER<VECTOR<float,2> >;
template class TEMPERATURE_CONTAINER<VECTOR<float,3> >;
template class TEMPERATURE_CONTAINER<VECTOR<double,1> >;
template class TEMPERATURE_CONTAINER<VECTOR<double,2> >;
template class TEMPERATURE_CONTAINER<VECTOR<double,3> >;
}
