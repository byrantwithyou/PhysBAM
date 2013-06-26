//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Sergey Koltakov, Neil Molino, Igor Neverov, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <Incompressible/Level_Sets/LEVELSET_COLLIDABLE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVELSET_COLLIDABLE<TV>::
LEVELSET_COLLIDABLE(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input,GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collision_body_list_input,const int number_of_ghost_cells_input,
    const bool set_secondary_interpolation)
    :LEVELSET<TV>(grid_input,phi_input,number_of_ghost_cells_input),collision_body_list(&collision_body_list_input),clamp_phi_with_collision_bodies(true),
    collision_aware_interpolation_plus(0),collision_aware_interpolation_minus(0),collision_unaware_interpolation(0),collidable_phi_replacement_value((T)1e-5)
{
    collision_aware_interpolation_plus=new T_LINEAR_INTERPOLATION_SCALAR;
    collision_aware_interpolation_minus=new LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<TV>,T>(*collision_body_list,&valid_mask_current,collidable_phi_replacement_value);
    if(set_secondary_interpolation) secondary_interpolation=collision_aware_interpolation_minus;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_COLLIDABLE<TV>::
~LEVELSET_COLLIDABLE()
{
    assert(!collision_unaware_interpolation);
    delete collision_aware_interpolation_plus;
    delete collision_aware_interpolation_minus;
}
//#####################################################################
// Function Collision_Aware_Phi
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_COLLIDABLE<TV>::
Collision_Aware_Phi(const TV& location) const
{
    assert(collision_body_list);
    return LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<TV>,T>(*collision_body_list,&valid_mask_current,collidable_phi_replacement_value).Clamped_To_Array(grid,phi,location);
}
namespace PhysBAM{
template class LEVELSET_COLLIDABLE<VECTOR<float,1> >;
template class LEVELSET_COLLIDABLE<VECTOR<float,2> >;
template class LEVELSET_COLLIDABLE<VECTOR<float,3> >;
template class LEVELSET_COLLIDABLE<VECTOR<double,1> >;
template class LEVELSET_COLLIDABLE<VECTOR<double,2> >;
template class LEVELSET_COLLIDABLE<VECTOR<double,3> >;
}
