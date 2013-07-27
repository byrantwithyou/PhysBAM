//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_WRAPPER_COLLIDABLE_CELL
//#####################################################################
#ifndef __ADVECTION_WRAPPER_COLLIDABLE_CELL__
#define __ADVECTION_WRAPPER_COLLIDABLE_CELL__

#include <Tools/Advection/ADVECTION.h>
namespace PhysBAM{
template<class TV> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

// Assumes NESTED_ADVECTION is of type ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_ 
template<class TV,class T2,class T_NESTED_LOOKUP,class T_NESTED_ADVECTION,class T_FACE_LOOKUP_COLLIDABLE>
class ADVECTION_WRAPPER_COLLIDABLE_CELL:public ADVECTION<TV,T2>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef ARRAY<bool,FACE_INDEX<TV::m> > T_FACE_ARRAYS_BOOL;
    typedef typename T_FACE_LOOKUP_COLLIDABLE::template REBIND_NESTED_LOOKUP<T_NESTED_LOOKUP>::TYPE T_FACE_LOOKUP_COLLIDABLE_NESTED_LOOKUP;
public:

    T_NESTED_ADVECTION& nested_advection;
    const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list;
    const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask;

    ADVECTION_WRAPPER_COLLIDABLE_CELL(T_NESTED_ADVECTION& nested_advection_input,const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list_input,const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input)
        :nested_advection(nested_advection_input),body_list(body_list_input),face_velocities_valid_mask(face_velocities_valid_mask_input)
    {}

    void Update_Advection_Equation_Cell_Lookup(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const T_NESTED_LOOKUP& face_velocities,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max)
    {T_FACE_LOOKUP_COLLIDABLE_NESTED_LOOKUP V_lookup(face_velocities,body_list,&face_velocities_valid_mask);
    nested_advection.Update_Advection_Equation_Cell_Lookup(grid,Z,Z_ghost,V_lookup,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);}

//#####################################################################
};
}
#endif
