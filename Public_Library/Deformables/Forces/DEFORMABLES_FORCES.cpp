//#####################################################################
// Copyright 2007-2008, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Vectors/VECTOR.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLES_FORCES<TV>::
DEFORMABLES_FORCES(DEFORMABLE_PARTICLES<TV>& particles)
    :particles(particles),cfl_number((T)1),allow_external_cfl_number(true),cfl_initialized(false),
    unique_id(Get_Unique_Id()),compute_half_forces(false)
{
    Use_Rest_State_For_Strain_Rate(false);
    Limit_Time_Step_By_Strain_Rate();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLES_FORCES<TV>::
~DEFORMABLES_FORCES()
{}
//#####################################################################
// Function Use_Rest_State_For_Strain_Rate
//#####################################################################
template<class TV> void DEFORMABLES_FORCES<TV>::
Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input)
{
    use_rest_state_for_strain_rate=use_rest_state_for_strain_rate_input;
}
//#####################################################################
// Function Limit_Time_Step_By_Strain_Rate
//#####################################################################
template<class TV> void DEFORMABLES_FORCES<TV>::
Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input,const T max_strain_per_time_step_input)
{
    limit_time_step_by_strain_rate=limit_time_step_by_strain_rate_input;
    assert(max_strain_per_time_step_input>0);max_strain_per_time_step=max_strain_per_time_step_input;
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int DEFORMABLES_FORCES<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void DEFORMABLES_FORCES<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void DEFORMABLES_FORCES<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void DEFORMABLES_FORCES<TV>::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void DEFORMABLES_FORCES<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void DEFORMABLES_FORCES<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLES_FORCES<TV>::
Potential_Energy(const T time) const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return 0;
}
//#####################################################################
// Function Add_Force_Data
//#####################################################################
template<class TV> void DEFORMABLES_FORCES<TV>::
Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name) const
{
};
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void DEFORMABLES_FORCES<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
}
//#####################################################################
// Function Test_Diff
//#####################################################################
// template<class TV> void SURFACE_TENSION_FORCE_3D<TV>::
// Test_Diff(const T time) 
// {
//     PHYSBAM_ASSERT(sizeof(T)==sizeof(double));
//     RANDOM_NUMBERS<T> random;
//     T e=(T)1e-5;
//     ARRAY<TV> dX(particles.X.m);
//     random.Fill_Uniform(dX,-e,e);
//     ARRAY<TV> X2a(particles.X+dX);
//     ARRAY_VIEW<TV> X1(X2a);
//     ARRAY<TV> F0(particles.X.m);
//     ARRAY<TV> F1(particles.X.m);
//     ARRAY<TV> G0(particles.X.m);
//     ARRAY<TV> G1(particles.X.m);
//     Update_Position_Based_State(time,true,false);
//     T psi0=Potential_Energy(time);
//     Add_Velocity_Independent_Forces(F0,time);
//     Add_Implicit_Velocity_Independent_Forces(dX,G0,time);
//     particles.X.Exchange(X1);
//     Update_Position_Based_State(time,true,false);
//     T psi1=Potential_Energy(time);
//     Add_Velocity_Independent_Forces(F1,time);
//     Add_Implicit_Velocity_Independent_Forces(dX,G1,time);
//     particles.X.Exchange(X1);
//     Update_Position_Based_State(time,true,false);
        
//     ARRAY<TV> F0pF1(particles.X.m);
//     ARRAY<TV> F1mF0(particles.X.m);
//     ARRAY<TV> G0pG1o2(particles.X.m);
//     for(int k=0;k<F0.m;k++){
//         F0pF1(k)=F0(k)+F1(k);
//         G0pG1o2(k)=(G0(k)+G1(k))/2;
//         F1mF0(k)=F1(k)-F0(k);}
//     T df=(F0pF1.Dot(dX)/2);
//     T dpsi=(psi1-psi0);
//     T error=(dpsi+df)/e;
//     LOG::cout<<"Energy Diff Test: "<<psi0<<" "<<psi1<<" "<<error<<std::endl;        
//     T MF=sqrt(F1mF0.Magnitude_Squared());
//     T MG=sqrt(G0pG1o2.Magnitude_Squared());
//     T ferror=sqrt((F1mF0-G0pG1o2).Magnitude_Squared());
//     LOG::cout<<"Force Diff Test: "<<MF<<" "<<MG<<" "<<ferror<<" "<<ferror/max((T)1e-30,MF,MG)<<std::endl;
// }
//#####################################################################
namespace PhysBAM{
#define INSTANTIATION_HELPER(T,d)                       \
    template class DEFORMABLES_FORCES<VECTOR<T,d> >;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
}
