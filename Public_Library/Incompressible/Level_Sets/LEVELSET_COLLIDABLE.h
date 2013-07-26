//#####################################################################
// Copyright 2002-2006, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_COLLIDABLE
//#####################################################################
#ifndef __LEVELSET_COLLIDABLE__
#define __LEVELSET_COLLIDABLE__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Tools/Math_Tools/constants.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <Tools/Vectors/SCALAR_POLICY.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Incompressible/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY_UNIFORM.h>
#include <cassert>
#include <cfloat>
namespace PhysBAM{

template<class TV> struct INTERPOLATION_POLICY;
template<class TV> class LEVELSET_CALLBACKS; // TODO: invalid dependency
template<class TV> struct BOUNDARY_POLICY;
template<class TV> struct GRID_ARRAYS_POLICY;
template<class TV> class GRID;
template<class TV,class T2> class BOUNDARY;
template<class TV> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

template<class TV>
class LEVELSET_COLLIDABLE:public LEVELSET<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef INTERPOLATION_UNIFORM<TV,T> T_INTERPOLATION_SCALAR;
public:
    using LEVELSET<TV>::interpolation;using LEVELSET<TV>::valid_mask_current;using LEVELSET<TV>::grid;using LEVELSET<TV>::phi;
    using LEVELSET<TV>::secondary_interpolation;

    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>* collision_body_list;
    bool collision_aware_signed_distance,clamp_phi_with_collision_bodies;
    T_INTERPOLATION_SCALAR *collision_aware_interpolation_plus,*collision_aware_interpolation_minus,*collision_unaware_interpolation;
    T collidable_phi_replacement_value;

    LEVELSET_COLLIDABLE(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input,GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_body_list_input,const int number_of_ghost_cells_input,
        const bool set_secondary_interpolation=false);
    ~LEVELSET_COLLIDABLE();

    void Enable_Collision_Aware_Interpolation(const int sign)
    {PHYSBAM_ASSERT(!collision_unaware_interpolation);collision_unaware_interpolation=interpolation;
    interpolation=sign==1?collision_aware_interpolation_plus:collision_aware_interpolation_minus;}

    void Disable_Collision_Aware_Interpolation()
    {interpolation=collision_unaware_interpolation;collision_unaware_interpolation=0;}

    T Collision_Aware_Phi(const TV& location) const;
//#####################################################################
};
}
#endif
