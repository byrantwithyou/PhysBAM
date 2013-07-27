//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TEMPERATURE_CONTAINER
//#####################################################################
#ifndef __TEMPERATURE_CONTAINER__
#define __TEMPERATURE_CONTAINER__

#include <Incompressible/Advection_Collidable/ADVECTION_COLLIDABLE_FORWARD.h>
#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_UNIFORM.h>
#include <Incompressible/Grid_Based_Fields/GRID_AND_ARRAY_CONTAINER.h>
namespace PhysBAM{

template<class TV> struct GRID_ARRAYS_POLICY;

template<class TV>
class TEMPERATURE_CONTAINER:public GRID_AND_ARRAY_CONTAINER<TV,typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename ADVECTION_COLLIDABLE_POLICY<TV>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;
    typedef FACE_LOOKUP_COLLIDABLE_UNIFORM<TV> T_FACE_LOOKUP_COLLIDABLE;
    typedef typename REBIND<ARRAY<T,FACE_INDEX<TV::m> >,bool>::TYPE T_FACE_ARRAYS_BOOL;
public:
    typedef GRID_AND_ARRAY_CONTAINER<TV,T> BASE;
    using BASE::grid;using BASE::array;using BASE::boundary_default;using BASE::boundary;using BASE::Set_To_Constant_Value;using BASE::Set_Custom_Advection;

    ARRAY<T,TV_INT>& temperature;
    T ambient_temperature;
    T cooling_constant;
    T hot_point;

private:
    ARRAY<bool,TV_INT> valid_mask_current;
    ARRAY<bool,TV_INT> valid_mask_next;
public:
    T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL* nested_semi_lagrangian_collidable;
    ADVECTION_WRAPPER_COLLIDABLE_CELL<TV,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL,T_FACE_LOOKUP_COLLIDABLE>* semi_lagrangian_collidable;

    TEMPERATURE_CONTAINER(GRID<TV>& grid_input);
    ~TEMPERATURE_CONTAINER();

    void Initialize_Array(const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true) PHYSBAM_OVERRIDE
    {GRID_AND_ARRAY_CONTAINER<TV,T>::Initialize_Array(ghost_cells,initialize_new_elements,copy_existing_elements);
    if(semi_lagrangian_collidable){valid_mask_current.Resize(grid.Cell_Indices(3),true,true,true);valid_mask_next.Resize(grid.Cell_Indices(3),false);}}
  
    void Set_Ambient_Temperature(const T temperature_input=298)
    {ambient_temperature=temperature_input;boundary->Set_Fixed_Boundary(true,ambient_temperature);}

    void Set_To_Ambient_Temperature()
    {Set_To_Constant_Value(ambient_temperature);}

    void Set_Cooling_Constant(const T cooling_constant_input=3000)
    {cooling_constant=cooling_constant_input;}

    void Set_Hot_Point(const T hot_point_input=3000)
    {hot_point=hot_point_input;}

//#####################################################################
    void Euler_Step(const T dt,const T time,const int number_of_ghost_cells) PHYSBAM_OVERRIDE;
    void Apply_Cooling(const T dt,const T time);
    void Apply_Individual_Cooling(T& temperature,const T constant);
    void Use_Semi_Lagrangian_Collidable_Advection(const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list,const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input);
//#####################################################################
};
}
#endif
