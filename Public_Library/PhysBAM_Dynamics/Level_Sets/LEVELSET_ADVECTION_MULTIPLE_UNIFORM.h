//#####################################################################
// Copyright 2009, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_ADVECTION_MULTIPLE_UNIFORM
//##################################################################### 
#ifndef __LEVELSET_ADVECTION_MULTIPLE_UNIFORM__
#define __LEVELSET_ADVECTION_MULTIPLE_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_MULTIPLE.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_UNIFORM.h>

namespace PhysBAM {
    
template<class T_GRID>
class LEVELSET_ADVECTION_MULTIPLE_UNIFORM:public LEVELSET_ADVECTION_MULTIPLE<T_GRID>
{
    typedef LEVELSET_ADVECTION_MULTIPLE<T_GRID> BASE;
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef typename ADVECTION_COLLIDABLE_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID> T_FACE_LOOKUP_COLLIDABLE;
    typedef typename REBIND<ARRAY<T,FACE_INDEX<TV::m> >,bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<T_GRID> T_FACE_LOOKUP_COLLIDABLE_SLIP;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET_MULTIPLE T_LEVELSET_MULTIPLE;
public:
    
    LEVELSET_ADVECTION_MULTIPLE_UNIFORM(T_LEVELSET_MULTIPLE* _levelsets):BASE(_levelsets){};
};

}
#endif
