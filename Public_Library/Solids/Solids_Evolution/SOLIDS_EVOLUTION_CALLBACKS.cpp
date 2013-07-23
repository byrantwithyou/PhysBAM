//#####################################################################
// Copyright 2004-2008, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <Solids/Solids_Evolution/SOLIDS_EVOLUTION_CALLBACKS.h>
#include <cfloat>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_EVOLUTION_CALLBACKS<TV>::
~SOLIDS_EVOLUTION_CALLBACKS()
{
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION_CALLBACKS<TV>::
Update_Solids_Parameters(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Self_Collisions_Begin_Callback
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION_CALLBACKS<TV>::
Self_Collisions_Begin_Callback(const T time,const int substep)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Preprocess_Solids_Substep
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION_CALLBACKS<TV>::
Preprocess_Solids_Substep(const T time,const int substep)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Postprocess_Solids_Substep
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION_CALLBACKS<TV>::
Postprocess_Solids_Substep(const T time,const int substep)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Align_Deformable_Bodies_With_Rigid_Bodies
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION_CALLBACKS<TV>::
Align_Deformable_Bodies_With_Rigid_Bodies()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Apply_Constraints
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION_CALLBACKS<TV>::
Apply_Constraints(const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Constraints_CFL
//#####################################################################
template<class TV> typename TV::SCALAR SOLIDS_EVOLUTION_CALLBACKS<TV>::
Constraints_CFL()
{
    return FLT_MAX;
}
//#####################################################################
// Function Limit_Solids_Dt
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION_CALLBACKS<TV>::
Limit_Solids_Dt(T& dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Pre_Advance_Cluster_Fracture
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION_CALLBACKS<TV>::
Pre_Advance_Cluster_Fracture(const T& dt, const T& time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Post_Advance_Cluster_Fracture
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION_CALLBACKS<TV>::
Post_Advance_Cluster_Fracture(const T& dt, const T& time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Filter_Velocities
//#####################################################################
template<class TV> void SOLIDS_EVOLUTION_CALLBACKS<TV>::
Filter_Velocities(const T dt,const T time,const bool velocity_update)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Solid_Source_Velocities
//#####################################################################
template<class TV> bool SOLIDS_EVOLUTION_CALLBACKS<TV>::
Get_Solid_Source_Velocities(ARRAY<int>& deformable_simplices,ARRAY<T>& deformable_simplex_forces,ARRAY<PAIR<int,int> >& rigid_simplices,ARRAY<T>& rigid_simplex_forces,TV& orientation,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return false;
}
namespace PhysBAM{
template class SOLIDS_EVOLUTION_CALLBACKS<VECTOR<double,1> >;
template class SOLIDS_EVOLUTION_CALLBACKS<VECTOR<double,2> >;
template class SOLIDS_EVOLUTION_CALLBACKS<VECTOR<double,3> >;
template class SOLIDS_EVOLUTION_CALLBACKS<VECTOR<float,1> >;
template class SOLIDS_EVOLUTION_CALLBACKS<VECTOR<float,2> >;
template class SOLIDS_EVOLUTION_CALLBACKS<VECTOR<float,3> >;
}
