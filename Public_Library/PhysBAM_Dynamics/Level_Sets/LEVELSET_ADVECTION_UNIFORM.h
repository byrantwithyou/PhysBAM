//#####################################################################
// Copyright 2009, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_ADVECTION_UNIFORM
//##################################################################### 
#ifndef __LEVELSET_ADVECTION_UNIFORM__
#define __LEVELSET_ADVECTION_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION.h>

namespace PhysBAM {
    
template<class T_GRID>
class LEVELSET_ADVECTION_UNIFORM:public LEVELSET_ADVECTION<T_GRID>
{
public:
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
private:
    typedef LEVELSET_ADVECTION<T_GRID> BASE;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef typename ADVECTION_COLLIDABLE_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE T_FACE_LOOKUP_COLLIDABLE;
    typedef typename REBIND<ARRAY<T,FACE_INDEX<TV::m> >,bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE_SLIP T_FACE_LOOKUP_COLLIDABLE_SLIP;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
public:
    using BASE::advection;
    using BASE::levelset;
    using BASE::reinitialization_cfl;using BASE::reinitialization_runge_kutta_order;using BASE::reinitialization_spatial_order;using BASE::Set_Custom_Advection;

    ADVECTION_MACCORMACK_UNIFORM<T_GRID,T,ADVECTION<T_GRID,T> >* advection_maccormack;

    LEVELSET_ADVECTION_UNIFORM(T_LEVELSET* _levelset)
        :BASE(_levelset),advection_maccormack(0)
    {}

    void Use_Maccormack_Advection(const ARRAY<bool,TV_INT>& cell_mask);
    T Approximate_Negative_Material(const T interface_thickness=3,const T time=0) const;
    T Approximate_Positive_Material(const T interface_thickness=3,const T time=0) const;
    void Euler_Step_Of_Reinitialization(const ARRAY<T,TV_INT>& sign_phi,const T dt,const T time);
    void Reinitialize(const int time_steps,const T time);
    void Euler_Step(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocity,const T dt,const T time,const int number_of_ghost_cells);
};
}
#endif
