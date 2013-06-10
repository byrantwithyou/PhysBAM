//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_1D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_FORCE_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <climits>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_FORCE_COLLECTION<TV>::
RIGID_FORCE_COLLECTION(RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
    :rigid_body_collection(rigid_body_collection),print_energy(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_FORCE_COLLECTION<TV>::
~RIGID_FORCE_COLLECTION()
{
    rigids_forces.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_FORCE_COLLECTION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const
{
    for(int k=0;k<rigids_forces.m;k++)
        if(rigids_forces(k)->use_velocity_independent_forces) rigids_forces(k)->Add_Velocity_Independent_Forces(rigid_F_full,time);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
// can depend on position too
template<class TV> void RIGID_FORCE_COLLECTION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const
{
    for(int k=0;k<rigids_forces.m;k++)
        if(rigids_forces(k)->use_velocity_dependent_forces) rigids_forces(k)->Add_Velocity_Dependent_Forces(rigid_V_full,rigid_F_full,time);
}
//#####################################################################
// Function Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_FORCE_COLLECTION<TV>::
Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T scale,const T time) const
{
    assert(rigid_F_full.Size()==rigid_body_collection.rigid_body_particles.Size());
    rigid_F_full.Subset(rigid_body_collection.dynamic_rigid_body_particles).Fill(TWIST<TV>()); // note we zero here because we will scale the forces below
    bool added=false;
    for(int k=0;k<rigids_forces.m;k++) if(rigids_forces(k)->use_implicit_velocity_independent_forces){
        rigids_forces(k)->Add_Implicit_Velocity_Independent_Forces(rigid_V_full,rigid_F_full,time);added=true;}
    if(added) rigid_F_full.Subset(rigid_body_collection.simulated_rigid_body_particles)*=scale;
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void RIGID_FORCE_COLLECTION<TV>::
Update_Position_Based_State(const T time)
{
    for(int k=0;k<rigids_forces.m;k++)
        if(rigids_forces(k)->use_position_based_state) rigids_forces(k)->Update_Position_Based_State(time);
}
//#####################################################################
// Function Compute_Energy
//#####################################################################
template<class TV> void RIGID_FORCE_COLLECTION<TV>::
Compute_Energy(const T time,T& kinetic_energy,T& potential_energy) const
{
    potential_energy=0;
    kinetic_energy=0;
    for(int i=0;i<rigids_forces.m;i++) potential_energy+=rigids_forces(i)->Potential_Energy(time);
    for(int i=0;i<rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=rigid_body_collection.dynamic_rigid_body_particles(i);
        kinetic_energy+=rigid_body_collection.Rigid_Body(p).Kinetic_Energy();}
}
//#####################################################################
// Function Print_Energy
//#####################################################################
template<class TV> void RIGID_FORCE_COLLECTION<TV>::
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
template<class TV> typename TV::SCALAR RIGID_FORCE_COLLECTION<TV>::
CFL_Rigid(const RIGID_BODY_EVOLUTION_PARAMETERS<TV>& rigid_body_evolution_parameters,const bool verbose_dt)
{
    static T static_min_bounding_box_width=FLT_MAX;
    T min_bounding_box_width=FLT_MAX;
    for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++) if(rigid_body_collection.Is_Active(i)){
            const RANGE<TV>& box=rigid_body_collection.Rigid_Body(i).Object_Space_Bounding_Box();
            TV edge_lengths=box.Edge_Lengths();min_bounding_box_width=min(min_bounding_box_width,edge_lengths.Min());}
    if(min_bounding_box_width!=static_min_bounding_box_width){
        static_min_bounding_box_width=min_bounding_box_width;
        LOG::Stat("minimum rigid body bounding box width",min_bounding_box_width);}

    T max_distance_per_time_step=rigid_body_evolution_parameters.max_rigid_body_linear_movement_fraction_per_time_step*min_bounding_box_width;
    T dt=FLT_MAX;
    bool no_active_bodies=true;
    for(int p=0;p<rigid_body_collection.rigid_body_particles.Size();p++) if(rigid_body_collection.Is_Active(p)){
        dt=min(dt,rigid_body_collection.Rigid_Body(p).CFL(max_distance_per_time_step,rigid_body_evolution_parameters.max_rigid_body_rotation_per_time_step,verbose_dt));
        no_active_bodies=false;}
    if(no_active_bodies) return FLT_MAX; // don't apply rigid dt bounds if there aren't any active rigid bodies
    dt=Robust_Multiply(rigid_body_evolution_parameters.rigid_cfl,dt);
    T dt_clamped=clamp(dt,rigid_body_evolution_parameters.rigid_minimum_dt,rigid_body_evolution_parameters.rigid_maximum_dt);
    if(dt_clamped>dt && verbose_dt) LOG::cout<<"Warning: taking larger time step ("<<dt_clamped<<") than CFL dt ("<<dt<<")"<<std::endl;
    return dt_clamped;
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int RIGID_FORCE_COLLECTION<TV>::
Add_Force(RIGIDS_FORCES<TV>* force)
{
    rigids_forces.Append(force);
    force->Set_CFL_Number((T).5);
    return rigids_forces.m;
}
template class RIGID_FORCE_COLLECTION<VECTOR<float,1> >;
template class RIGID_FORCE_COLLECTION<VECTOR<float,2> >;
template class RIGID_FORCE_COLLECTION<VECTOR<float,3> >;
template class RIGID_FORCE_COLLECTION<VECTOR<double,1> >;
template class RIGID_FORCE_COLLECTION<VECTOR<double,2> >;
template class RIGID_FORCE_COLLECTION<VECTOR<double,3> >;
}
