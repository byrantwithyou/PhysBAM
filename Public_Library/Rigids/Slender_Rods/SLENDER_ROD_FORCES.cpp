//#####################################################################
// Copyright 2007-2008, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Vectors/VECTOR.h>
#include <Rigids/Slender_Rods/SLENDER_ROD_FORCES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SLENDER_ROD_FORCES<TV>::
SLENDER_ROD_FORCES(RIGID_SLENDER_ROD_COLLECTION<TV>& rigid_slender_rod_collection)
    :rigid_slender_rod_collection(rigid_slender_rod_collection),cfl_number((T)1),allow_external_cfl_number(true),cfl_initialized(false),
    unique_id(Get_Unique_Id()),compute_half_forces(false)
{
    Use_Rest_State_For_Strain_Rate(false);
    Limit_Time_Step_By_Strain_Rate();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SLENDER_ROD_FORCES<TV>::
~SLENDER_ROD_FORCES()
{}
//#####################################################################
// Function Use_Rest_State_For_Strain_Rate
//#####################################################################
template<class TV> void SLENDER_ROD_FORCES<TV>::
Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input)
{
    use_rest_state_for_strain_rate=use_rest_state_for_strain_rate_input;
}
//#####################################################################
// Function Limit_Time_Step_By_Strain_Rate
//#####################################################################
template<class TV> void SLENDER_ROD_FORCES<TV>::
Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input,const T max_strain_per_time_step_input)
{
    limit_time_step_by_strain_rate=limit_time_step_by_strain_rate_input;
    assert(max_strain_per_time_step_input>0);max_strain_per_time_step=max_strain_per_time_step_input;
}
template<class TV> int SLENDER_ROD_FORCES<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void SLENDER_ROD_FORCES<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void SLENDER_ROD_FORCES<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> void SLENDER_ROD_FORCES<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time,bool transpose) const
{
}
template<class TV> void SLENDER_ROD_FORCES<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class TV> typename TV::SCALAR SLENDER_ROD_FORCES<TV>::
Potential_Energy(const T time) const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return 0;
}
//#####################################################################
// Function Update_Force_Elements
//#####################################################################
template<int d,class T_ARRAY> void
Update_Force_Elements(ARRAY<int>& force_elements,
    const ARRAY_BASE<VECTOR<int,d>,T_ARRAY>& elements,
    const ARRAY<bool>& particle_is_simulated)
{
    const T_ARRAY& e=elements.Derived();
    force_elements.Remove_All();
    for(int i=0;i<e.Size();i++)
        if(particle_is_simulated.Subset(e(i)).Contains(true))
            force_elements.Append(i);
}
//#####################################################################
// Function Update_Force_Particles
//#####################################################################
template<class T_ARRAY> void
Update_Force_Particles(ARRAY<int>& force_particles,
    const ARRAY_BASE<int,T_ARRAY>& particles,
    const ARRAY<bool>& particle_is_simulated,bool check_dups)
{
    const T_ARRAY& e=particles.Derived();
    force_particles.Remove_All();
    for(int i=0;i<e.Size();i++)
        if(particle_is_simulated(e(i)))
            force_particles.Append(e(i));
    if(check_dups) Prune_Duplicates(force_particles);
}
//#####################################################################
template class SLENDER_ROD_FORCES<VECTOR<float,1> >;
template class SLENDER_ROD_FORCES<VECTOR<float,2> >;
template class SLENDER_ROD_FORCES<VECTOR<float,3> >;
template class SLENDER_ROD_FORCES<VECTOR<double,1> >;
template class SLENDER_ROD_FORCES<VECTOR<double,2> >;
template class SLENDER_ROD_FORCES<VECTOR<double,3> >;
}
