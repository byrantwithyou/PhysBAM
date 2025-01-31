//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_WRAPPER_COLLIDABLE_FACE
//#####################################################################
#ifndef __ADVECTION_WRAPPER_COLLIDABLE_FACE__
#define __ADVECTION_WRAPPER_COLLIDABLE_FACE__

#include <Grid_PDE/Advection/ADVECTION.h>

namespace PhysBAM{

template<class TV> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

// Assumes NESTED_ADVECTION is of type ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE
template<class TV,class T2,class T_NESTED_LOOKUP,class T_NESTED_ADVECTION,class T_FACE_LOOKUP_COLLIDABLE>
class ADVECTION_WRAPPER_COLLIDABLE_FACE:public ADVECTION<TV,T2,T_NESTED_LOOKUP>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename T_FACE_LOOKUP_COLLIDABLE::template REBIND_NESTED_LOOKUP<T_NESTED_LOOKUP>::TYPE T_FACE_LOOKUP_COLLIDABLE_NESTED_LOOKUP;
public:
    T_NESTED_ADVECTION& nested_advection;
    const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list;
    const ARRAY<bool,FACE_INDEX<TV::m> >& face_velocities_valid_mask;

    ADVECTION_WRAPPER_COLLIDABLE_FACE(T_NESTED_ADVECTION& nested_advection_input,GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list_input,const ARRAY<bool,FACE_INDEX<TV::m> >& face_velocities_valid_mask_input)
        :nested_advection(nested_advection_input),body_list(body_list_input),face_velocities_valid_mask(face_velocities_valid_mask_input)
    {}

    void Update_Advection_Equation_Face_Lookup(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& Z,const T_NESTED_LOOKUP& Z_ghost,
        const T_NESTED_LOOKUP& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
        const T_NESTED_LOOKUP* Z_min_ghost,const T_NESTED_LOOKUP* Z_max_ghost,ARRAY<T,FACE_INDEX<TV::m> >* Z_min,ARRAY<T,FACE_INDEX<TV::m> >* Z_max)
    {T_FACE_LOOKUP_COLLIDABLE_NESTED_LOOKUP Z_ghost_lookup(Z_ghost,body_list,&face_velocities_valid_mask);
    T_FACE_LOOKUP_COLLIDABLE_NESTED_LOOKUP V_lookup(face_velocities,body_list,&face_velocities_valid_mask);
    if(Z_min_ghost && Z_max_ghost){
        T_FACE_LOOKUP_COLLIDABLE_NESTED_LOOKUP Z_min_ghost_lookup(*Z_min_ghost,body_list,&face_velocities_valid_mask),Z_max_ghost_lookup(*Z_max_ghost,body_list,&face_velocities_valid_mask);
        nested_advection.Update_Advection_Equation_Face_Lookup(grid,Z,Z_ghost_lookup,V_lookup,boundary,dt,time,&Z_min_ghost_lookup,&Z_max_ghost_lookup,Z_min,Z_max);}
    else nested_advection.Update_Advection_Equation_Face_Lookup(grid,Z,Z_ghost_lookup,V_lookup,boundary,dt,time,0,0,Z_min,Z_max);}

//#####################################################################
};
}
#endif
