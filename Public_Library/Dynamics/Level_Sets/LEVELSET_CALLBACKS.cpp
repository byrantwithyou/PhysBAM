//#####################################################################
// Copyright 2003-2005, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_CALLBACKS<TV>::
~LEVELSET_CALLBACKS()
{
}
//#####################################################################
// Function Get_Levelset_Velocity
//#####################################################################
template<class TV> void LEVELSET_CALLBACKS<TV>::
Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET<TV>& levelset,ARRAY<T,FACE_INDEX<TV::m> >& face_velocity,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Get_Levelset_Velocity
//#####################################################################
template<class TV> void LEVELSET_CALLBACKS<TV>::
Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET_MULTIPLE<TV>& levelset,ARRAY<T,FACE_INDEX<TV::m> >& face_velocity,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Adjust_Particle_For_Domain_Boundaries
//#####################################################################
template<class TV> void LEVELSET_CALLBACKS<TV>::
Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,
    const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
template<class TV> bool LEVELSET_CALLBACKS<TV>::
Adjust_Particle_For_Objects(TV& X,TV& V,const T r, const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return true;
}
//#####################################################################
namespace PhysBAM{
template class LEVELSET_CALLBACKS<VECTOR<double,1> >;
template class LEVELSET_CALLBACKS<VECTOR<double,2> >;
template class LEVELSET_CALLBACKS<VECTOR<double,3> >;
template class LEVELSET_CALLBACKS<VECTOR<float,1> >;
template class LEVELSET_CALLBACKS<VECTOR<float,2> >;
template class LEVELSET_CALLBACKS<VECTOR<float,3> >;
}
