//#####################################################################
// Copyright 2007-2008, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Vectors/VECTOR.h>
#include <Solids/Forces_And_Torques/SOLIDS_FORCES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_FORCES<TV>::
SOLIDS_FORCES(DEFORMABLE_PARTICLES<TV>& particles,RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
    :particles(particles),rigid_body_collection(rigid_body_collection),cfl_number((T)1),allow_external_cfl_number(true),
    cfl_initialized(false),unique_id(Get_Unique_Id()),compute_half_forces(false)
{
    Use_Rest_State_For_Strain_Rate(false);
    Limit_Time_Step_By_Strain_Rate();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_FORCES<TV>::
~SOLIDS_FORCES()
{}
//#####################################################################
// Function Use_Rest_State_For_Strain_Rate
//#####################################################################
template<class TV> void SOLIDS_FORCES<TV>::
Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input)
{
    use_rest_state_for_strain_rate=use_rest_state_for_strain_rate_input;
}
//#####################################################################
// Function Limit_Time_Step_By_Strain_Rate
//#####################################################################
template<class TV> void SOLIDS_FORCES<TV>::
Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input,const T max_strain_per_time_step_input)
{
    limit_time_step_by_strain_rate=limit_time_step_by_strain_rate_input;
    assert(max_strain_per_time_step_input>0);max_strain_per_time_step=max_strain_per_time_step_input;
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int SOLIDS_FORCES<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void SOLIDS_FORCES<TV>::
Add_Velocity_Dependent_Forces_First_Half(const GENERALIZED_VELOCITY<TV>& V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void SOLIDS_FORCES<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,GENERALIZED_VELOCITY<TV>& F,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SOLIDS_FORCES<TV>::
Add_Implicit_Velocity_Independent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time,bool transpose) const
{
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void SOLIDS_FORCES<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR SOLIDS_FORCES<TV>::
Potential_Energy(const T time) const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return 0;
}
//#####################################################################
// Function Add_Force_Data
//#####################################################################
template<class TV> void SOLIDS_FORCES<TV>::
Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name) const
{
};
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void SOLIDS_FORCES<TV>::
Update_Position_Based_State(const T time)
{
}
//#####################################################################
namespace PhysBAM{
#define INSTANTIATION_HELPER(T,d) \
    template class SOLIDS_FORCES<VECTOR<T,d> >;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
}
