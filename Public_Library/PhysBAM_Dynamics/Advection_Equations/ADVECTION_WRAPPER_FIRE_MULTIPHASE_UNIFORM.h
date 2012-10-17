//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_WRAPPER_FIRE_MULTIPHASE_UNIFORM
//#####################################################################
#ifndef __ADVECTION_WRAPPER_FIRE_MULTIPHASE_UNIFORM__
#define __ADVECTION_WRAPPER_FIRE_MULTIPHASE_UNIFORM__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_POLICY.h>
#include <PhysBAM_Dynamics/Interpolation/FIRE_INTERPOLATION_FORWARD.h>
namespace PhysBAM{

template<class T_GRID> class LEVELSET_MULTIPLE;

template<class T_GRID,class T2,class T_NESTED_ADVECTION>
class ADVECTION_WRAPPER_FIRE_MULTIPHASE_UNIFORM:public ADVECTION<T_GRID,T2>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef typename INCOMPRESSIBLE_POLICY<T_GRID>::PROJECTION T_PROJECTION;
public:
    T_NESTED_ADVECTION& nested_advection;
    const T_PROJECTION& projection;
    const LEVELSET_MULTIPLE<T_GRID>& levelset_multiple_n_plus_one;

    ADVECTION_WRAPPER_FIRE_MULTIPHASE_UNIFORM(T_NESTED_ADVECTION& nested_advection_input,const T_PROJECTION& projection_input,const LEVELSET_MULTIPLE<T_GRID>& levelset_multiple_n_plus_one_input)
        :nested_advection(nested_advection_input),projection(projection_input),levelset_multiple_n_plus_one(levelset_multiple_n_plus_one_input)
    {}

    void Update_Advection_Equation_Node(const T_GRID& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const ARRAY<TV,TV_INT>& V,BOUNDARY_UNIFORM<T_GRID,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0)
    {nested_advection.Update_Advection_Equation_Node(grid,Z,Z_ghost,V,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);}

    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const FACE_LOOKUP_UNIFORM<T_GRID>& face_velocities,BOUNDARY_UNIFORM<T_GRID,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max)
    {FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<T_GRID> V_lookup(face_velocities.V_face,projection,&levelset_multiple_n_plus_one);
    nested_advection.Update_Advection_Equation_Cell_Lookup(grid,Z,Z_ghost,V_lookup,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);}

    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const FACE_LOOKUP_UNIFORM<T_GRID>& Z_ghost,
        const FACE_LOOKUP_UNIFORM<T_GRID>& face_velocities,BOUNDARY_UNIFORM<T_GRID,T>& boundary,const T dt,const T time,
        const FACE_LOOKUP_UNIFORM<T_GRID>* Z_min_ghost,const FACE_LOOKUP_UNIFORM<T_GRID>* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
    {const LEVELSET_MULTIPLE<T_GRID>* levelset_multiple_n=projection.poisson_collidable->levelset_multiple; //assumes poisson's internal levelset is at time n
    FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<T_GRID> Z_ghost_lookup(Z_ghost.V_face,projection,levelset_multiple_n);
    FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<T_GRID> V_lookup(face_velocities.V_face,projection,&levelset_multiple_n_plus_one);
    if(Z_min_ghost && Z_max_ghost){
        FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<T_GRID> Z_min_ghost_lookup(Z_min_ghost->V_face,projection,&levelset_multiple_n_plus_one),
            Z_max_ghost_lookup(Z_max_ghost->V_face,projection,&levelset_multiple_n_plus_one);
        nested_advection.Update_Advection_Equation_Face_Lookup(grid,Z,Z_ghost_lookup,V_lookup,boundary,dt,time,&Z_min_ghost_lookup,&Z_max_ghost_lookup,Z_min,Z_max);}
    else nested_advection.Update_Advection_Equation_Face_Lookup(grid,Z,Z_ghost_lookup,V_lookup,boundary,dt,time,0,0,Z_min,Z_max);}

//#####################################################################
};
}
#endif
