//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/IDENTITY_ARRAY.h>
//#include <Core/Arrays/INDIRECT_ARRAY.h>
//#include <Core/Log/DEBUG_UTILITIES.h>
//#include <Core/Log/LOG.h>
//#include <Core/Read_Write/FILE_UTILITIES.h>
//#include <Core/Utilities/DEBUG_CAST.h>
//#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
//#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
//#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
//#include <Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
//#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
//#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Rigids/Slender_Rods/RIGID_SLENDER_ROD_COLLECTION.h>
#include <Rigids/Slender_Rods/RIGID_SLENDER_ROD_PARTICLES.h>
#include <Rigids/Slender_Rods/SLENDER_ROD_FORCES.h>
#include <climits>
//#include <Slender_Rod/Rigid_Bodies/RIGID_SLENDER_ROD.h>
//#include <Slender_Rod/Rigid_Bodies/RIGID_SLENDER_ROD_EVOLUTION_PARAMETERS.h>
//#include <Slender_Rod/Rigid_Slender_Rod_Clusters/RIGID_SLENDER_ROD_CLUSTER_BINDINGS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_SLENDER_ROD_COLLECTION<TV>::
RIGID_SLENDER_ROD_COLLECTION()
    :rigid_slender_rod_particles(*new RIGID_SLENDER_ROD_PARTICLES<TV>())
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_SLENDER_ROD_COLLECTION<TV>::
~RIGID_SLENDER_ROD_COLLECTION()
{
    delete &rigid_slender_rod_particles;
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void RIGID_SLENDER_ROD_COLLECTION<TV>::
Update_Simulated_Particles()
{
    simulated_rigid_rod_particles=IDENTITY_ARRAY<>(rigid_slender_rod_particles.number);
    dynamic_rigid_rod_particles=simulated_rigid_rod_particles;
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_SLENDER_ROD_COLLECTION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const
{
    for(int k=0;k<slender_rod_forces.m;k++)
        slender_rod_forces(k)->Add_Velocity_Independent_Forces(rigid_F_full,time);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
// can depend on position too
template<class TV> void RIGID_SLENDER_ROD_COLLECTION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const
{
    for(int k=0;k<slender_rod_forces.m;k++)
        slender_rod_forces(k)->Add_Velocity_Dependent_Forces(
            rigid_V_full,rigid_F_full,time);
}
//#####################################################################
// Function Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_SLENDER_ROD_COLLECTION<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,
    ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time,bool transpose) const
{
    assert(rigid_F_full.Size()==rigid_slender_rod_particles.Size());
    for(int k=0;k<slender_rod_forces.m;k++)
        slender_rod_forces(k)->Add_Implicit_Velocity_Independent_Forces(
            rigid_V_full,rigid_F_full,time,transpose);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void RIGID_SLENDER_ROD_COLLECTION<TV>::
Update_Position_Based_State(const T time)
{
    for(int k=0;k<slender_rod_forces.m;k++)
        slender_rod_forces(k)->Update_Position_Based_State(time);
}
//#####################################################################
// Function Compute_Energy
//#####################################################################
template<class TV> void RIGID_SLENDER_ROD_COLLECTION<TV>::
Compute_Energy(const T time,T& kinetic_energy,T& potential_energy) const
{
    potential_energy=0;
    kinetic_energy=0;
    for(int i=0;i<slender_rod_forces.m;i++)
        potential_energy+=slender_rod_forces(i)->Potential_Energy(time);
    for(int i=0;i<dynamic_rigid_rod_particles.m;i++){
//        int p=dynamic_rigid_rod_particles(i);
        PHYSBAM_NOT_IMPLEMENTED();
        //kinetic_energy+=Rigid_Slender_Rod(p).Kinetic_Energy();
    }
}
//#####################################################################
// Function Print_Energy
//#####################################################################
template<class TV> void RIGID_SLENDER_ROD_COLLECTION<TV>::
Print_Energy(const T time,const int step) const
{
    if(print_energy){
        T potential_energy=0,kinetic_energy=0;
        Compute_Energy(time,kinetic_energy,potential_energy);
        LOG::cout<<"total energy = "<<(potential_energy+kinetic_energy)<<"    (KE = "<<kinetic_energy<<"   PE = "<<potential_energy<<")  Step "<<step<<std::endl;}
}
//#####################################################################
// Function CFL_Rigid
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_SLENDER_ROD_COLLECTION<TV>::
CFL_Rigid_Rod(const bool verbose_dt)
{
    return FLT_MAX;
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int RIGID_SLENDER_ROD_COLLECTION<TV>::
Add_Force(SLENDER_ROD_FORCES<TV>* force)
{
    slender_rod_forces.Append(force);
    force->Set_CFL_Number((T).5);
    return slender_rod_forces.m;
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void RIGID_SLENDER_ROD_COLLECTION<TV>::
Read(const VIEWER_DIR& viewer_dir)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void RIGID_SLENDER_ROD_COLLECTION<TV>::
Write(const STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
template class RIGID_SLENDER_ROD_COLLECTION<VECTOR<float,1> >;
template class RIGID_SLENDER_ROD_COLLECTION<VECTOR<float,2> >;
template class RIGID_SLENDER_ROD_COLLECTION<VECTOR<float,3> >;
template class RIGID_SLENDER_ROD_COLLECTION<VECTOR<double,1> >;
template class RIGID_SLENDER_ROD_COLLECTION<VECTOR<double,2> >;
template class RIGID_SLENDER_ROD_COLLECTION<VECTOR<double,3> >;
}
