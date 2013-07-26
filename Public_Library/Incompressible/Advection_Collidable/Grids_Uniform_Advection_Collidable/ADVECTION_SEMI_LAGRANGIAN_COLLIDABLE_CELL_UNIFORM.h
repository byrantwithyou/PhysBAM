//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM
//#####################################################################
#ifndef __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM__
#define __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM__

#include <Tools/Advection/ADVECTION.h>
#include <Tools/Data_Structures/PAIR.h>
#include <Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_UNIFORM_FORWARD.h>
#include <Incompressible/Interpolation_Collidable/AVERAGING_COLLIDABLE_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
namespace PhysBAM{

template<class TV,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_UNIFORM<TV>
class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM:public ADVECTION<TV,T2,T_FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list;
    ARRAY<bool,TV_INT> &cell_valid_points_current,&cell_valid_points_next;
    T2 cell_crossover_replacement_value;
    bool extrapolate_to_revalidate_interpolation;
private:
    LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<TV,T2> linear_interpolation_collidable;
    LINEAR_INTERPOLATION_UNIFORM<TV,T2,typename T_FACE_LOOKUP::NESTED_LOOKUP> linear_interpolation;
    AVERAGING_UNIFORM<TV,typename T_FACE_LOOKUP::NESTED_LOOKUP> velocity_averaging;
    AVERAGING_COLLIDABLE_UNIFORM<TV,T_FACE_LOOKUP> velocity_averaging_collidable;
public:

    ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM(const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list_input,ARRAY<bool,TV_INT>& cell_valid_points_current_input,
        ARRAY<bool,TV_INT>& cell_valid_points_next_input,const T2& default_cell_replacement_value_input,const bool extrapolate_to_revalidate_interpolation_input);
    virtual ~ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM();

    T2 Compute_Revalidation_Value(const TV& from,const TV& to,const T2& current_invalid_value,const T2& default_value)
    {return default_value;}

//#####################################################################
    void Update_Advection_Equation_Cell_Lookup(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max);
    void Average_To_Invalidated_Cells(const GRID<TV>& grid,const T2 default_value,ARRAY<T2,TV_INT>& values);
//#####################################################################
};
}
#endif
