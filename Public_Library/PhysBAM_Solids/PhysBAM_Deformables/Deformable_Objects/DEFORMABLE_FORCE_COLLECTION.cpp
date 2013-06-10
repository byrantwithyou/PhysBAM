//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Levine, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_BODY_COLLECTION
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_FORCE_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLE_FORCE_COLLECTION<TV>::
DEFORMABLE_FORCE_COLLECTION(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection)
    :deformable_body_collection(deformable_body_collection),implicit_damping(true),print_energy(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLE_FORCE_COLLECTION<TV>::
~DEFORMABLE_FORCE_COLLECTION()
{
    deformables_forces.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Update_CFL
//#####################################################################
template<class TV> void DEFORMABLE_FORCE_COLLECTION<TV>::
Update_CFL()
{
    bool cfl_valid=true;
    if(deformables_forces.m){
        for(int i=0;i<deformables_forces.m;i++){if(!deformables_forces(i)->CFL_Valid()){cfl_valid=false;break;}}}
    else cfl_valid=false;
    if(!cfl_valid){
        frequency.Resize(deformable_body_collection.particles.Size(),false,false);
        INDIRECT_ARRAY<ARRAY<T_FREQUENCY_DEFORMABLE>,ARRAY<int>&> frequency_subset=frequency.Subset(deformable_body_collection.simulated_particles);
        frequency_subset.Fill(T_FREQUENCY_DEFORMABLE());

        for(int i=0;i<deformables_forces.m;i++){deformables_forces(i)->Initialize_CFL(frequency);deformables_forces(i)->Validate_CFL();}
        cfl_elastic=FLT_MAX;cfl_damping=FLT_MAX;
        for(int i=0;i<deformable_body_collection.simulated_particles.m;i++){int p=deformable_body_collection.simulated_particles(i);
            cfl_elastic=min(cfl_elastic,Robust_Inverse(sqrt(frequency(p).elastic_squared)));
            cfl_damping=min(cfl_damping,Robust_Inverse(frequency(p).damping));}}
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_FORCE_COLLECTION<TV>::
CFL(const bool verbose)
{
    T dt_elastic_and_damping=CFL_Elastic_And_Damping(),dt_strain_rate=CFL_Strain_Rate();
    if(verbose){
        LOG::cout<<"dt_elastic_and_damping = "<<dt_elastic_and_damping<<std::endl;
        LOG::cout<<"dt_strain_rate = "<<dt_strain_rate<<std::endl;
        LOG::cout<<"min = "<<min(dt_elastic_and_damping,dt_strain_rate)<<std::endl;}
    return min(dt_elastic_and_damping,dt_strain_rate);
}
//#####################################################################
// Function CFL_Elastic_And_Damping
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_FORCE_COLLECTION<TV>::
CFL_Elastic_And_Damping()
{
    T dt_elastic=CFL_Elastic();
    T dt_damping=FLT_MAX;if(!implicit_damping) dt_damping=CFL_Damping();
    T one_over_dt_full=1/dt_elastic+1/dt_damping;
    return Robust_Divide((T)1,one_over_dt_full);
}
//#####################################################################
// Function CFL_Elastic
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_FORCE_COLLECTION<TV>::
CFL_Elastic()
{
    Update_CFL();
    return cfl_elastic;
}
//#####################################################################
// Function CFL_Damping
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_FORCE_COLLECTION<TV>::
CFL_Damping()
{
    Update_CFL();
    return cfl_damping;
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_FORCE_COLLECTION<TV>::
CFL_Strain_Rate()
{
    T dt_strain=FLT_MAX;
    for(int k=0;k<deformables_forces.m;k++) if(deformables_forces(k)->limit_time_step_by_strain_rate) dt_strain=min(dt_strain,deformables_forces(k)->CFL_Strain_Rate()); // otherwise not included
    return dt_strain;
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void DEFORMABLE_FORCE_COLLECTION<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    for(int k=0;k<deformables_forces.m;k++) if(deformables_forces(k)->use_position_based_state) deformables_forces(k)->Update_Position_Based_State(time,is_position_update);
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_FORCE_COLLECTION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F_full,const T time) const
{
    assert(F_full.Size()==deformable_body_collection.particles.Size());
    for(int k=0;k<deformables_forces.m;k++) if(deformables_forces(k)->use_velocity_independent_forces) deformables_forces(k)->Add_Velocity_Independent_Forces(F_full,time);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
// can depend on position too
template<class TV> void DEFORMABLE_FORCE_COLLECTION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V_full,ARRAY_VIEW<TV> F_full,const T time) const
{
    assert(F_full.Size()==deformable_body_collection.particles.Size());
    for(int k=0;k<deformables_forces.m;k++) if(deformables_forces(k)->use_velocity_dependent_forces) deformables_forces(k)->Add_Velocity_Dependent_Forces(V_full,F_full,time);
}
//#####################################################################
// Function Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_FORCE_COLLECTION<TV>::
Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V_full,ARRAY_VIEW<TV> F_full,const T scale,const T time) const
{
    assert(V_full.Size()==deformable_body_collection.particles.Size() && F_full.Size()==deformable_body_collection.particles.Size());
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,ARRAY<int>&> F_subset=F_full.Subset(deformable_body_collection.dynamic_particles);
    F_subset.Fill(TV()); // note we zero here because we will scale the forces below
    bool added=false;
    for(int k=0;k<deformables_forces.m;k++) if(deformables_forces(k)->use_implicit_velocity_independent_forces){
        deformables_forces(k)->Add_Implicit_Velocity_Independent_Forces(V_full,F_full,time);added=true;}
    if(added) F_full.Subset(deformable_body_collection.simulated_particles)*=scale;
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int DEFORMABLE_FORCE_COLLECTION<TV>::
Add_Force(DEFORMABLES_FORCES<TV>* force)
{
    deformables_forces.Append(force);
    force->Set_CFL_Number(cfl_number);
    return deformables_forces.m;
}
template<class TV> void DEFORMABLE_FORCE_COLLECTION<TV>::
Set_CFL_Number(const T cfl_number_input)
{
    cfl_number=cfl_number_input;
    for(int i=0;i<deformables_forces.m;i++) deformables_forces(i)->Set_CFL_Number(cfl_number_input);
}
//#####################################################################
// Function Test_Energy
//#####################################################################
template<class TV> void DEFORMABLE_FORCE_COLLECTION<TV>::
Test_Energy(const T time)
{
    PHYSBAM_ASSERT(sizeof(T)==sizeof(double));
    RANDOM_NUMBERS<T> random;
    T e=(T)1e-5;
    ARRAY<TV> dX(deformable_body_collection.particles.X.m);
    random.Fill_Uniform(dX,-e,e);
    ARRAY<TV> X2a(deformable_body_collection.particles.X+dX);
    ARRAY_VIEW<TV> X2(X2a);
    for(int i=0;i<deformables_forces.m;i++){
        ARRAY<TV> F(deformable_body_collection.particles.X.m);
        deformables_forces(i)->Update_Position_Based_State(time,true);
        deformables_forces(i)->Add_Velocity_Independent_Forces(F,time);
        T PE1=deformables_forces(i)->Potential_Energy(time);
        deformable_body_collection.particles.X.Exchange(X2);
        deformables_forces(i)->Update_Position_Based_State(time,true);
        deformables_forces(i)->Add_Velocity_Independent_Forces(F,time);
        T PE2=deformables_forces(i)->Potential_Energy(time);
        deformable_body_collection.particles.X.Exchange(X2);
        deformables_forces(i)->Update_Position_Based_State(time,true);
        T W=F.Dot(dX)/2;
        T dPE=(PE1-PE2)/e,dW=W/e,rel=(dPE-dW)/max(abs(dW),(T)1e-20);
        LOG::cout<<"potential energy test d phi "<<dPE<<"  W "<<dW<<"   rel "<<rel<<"   "<<typeid(*deformables_forces(i)).name()<<std::endl;
    }
}
//#####################################################################
// Function Test_Force_Derivatives
//#####################################################################
template<class TV> void DEFORMABLE_FORCE_COLLECTION<TV>::
Test_Force_Derivatives(const T time)
{
    PHYSBAM_ASSERT(sizeof(T)==sizeof(double));
    RANDOM_NUMBERS<T> random;
    T e=(T)1e-5;
    ARRAY<TV> dX(deformable_body_collection.particles.X.m);
    random.Fill_Uniform(dX,-e,e);
    ARRAY<TV> X2a(deformable_body_collection.particles.X+dX);
    ARRAY_VIEW<TV> X2(X2a);
    for(int i=0;i<deformables_forces.m;i++){
        ARRAY<TV> F(deformable_body_collection.particles.X.m),G(deformable_body_collection.particles.X.m);
        deformables_forces(i)->Update_Position_Based_State(time,true);
        deformables_forces(i)->Add_Velocity_Independent_Forces(F,time);
        deformables_forces(i)->Add_Implicit_Velocity_Independent_Forces(dX,G,time);
        F*=-(T)1;
        deformable_body_collection.particles.X.Exchange(X2);
        deformables_forces(i)->Update_Position_Based_State(time,true);
        deformables_forces(i)->Add_Velocity_Independent_Forces(F,time);
        deformables_forces(i)->Add_Implicit_Velocity_Independent_Forces(dX,G,time);
        G*=(T).5;
        deformable_body_collection.particles.X.Exchange(X2);
        deformables_forces(i)->Update_Position_Based_State(time,true);
        T MF=sqrt(F.Magnitude_Squared());
        T MG=sqrt(G.Magnitude_Squared());
        T MD=sqrt((F-G).Magnitude_Squared());
        LOG::cout<<"force derivative error "<<MD<<" vs "<<MF<<"   "<<MG<<"   rel "<<MD/max((T)1e-30,MF,MG)<<"    "<<typeid(*deformables_forces(i)).name()<<std::endl;
    }
}
namespace PhysBAM{
template class DEFORMABLE_FORCE_COLLECTION<VECTOR<float,1> >;
template class DEFORMABLE_FORCE_COLLECTION<VECTOR<float,2> >;
template class DEFORMABLE_FORCE_COLLECTION<VECTOR<float,3> >;
template class DEFORMABLE_FORCE_COLLECTION<VECTOR<double,1> >;
template class DEFORMABLE_FORCE_COLLECTION<VECTOR<double,2> >;
template class DEFORMABLE_FORCE_COLLECTION<VECTOR<double,3> >;
}
