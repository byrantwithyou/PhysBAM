//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Levine, Igor Neverov, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_BODY_COLLECTION
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COLLISION_PENALTY_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_FORCE_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLID_BODY_COLLECTION<TV>::
SOLID_BODY_COLLECTION(EXAMPLE_FORCES_AND_VELOCITIES<TV>* example_forces_and_velocities_input)
    :collision_body_list(*new COLLISION_GEOMETRY_COLLECTION<TV>),
    deformable_body_collection(*new DEFORMABLE_BODY_COLLECTION<TV>(example_forces_and_velocities_input,collision_body_list)),
    rigid_body_collection(*new RIGID_BODY_COLLECTION<TV>(example_forces_and_velocities_input,&collision_body_list)),
    solid_force_collection(*new SOLID_FORCE_COLLECTION<TV>(deformable_body_collection.deformable_force_collection,rigid_body_collection.rigid_force_collection)),
    example_forces_and_velocities(example_forces_and_velocities_input),simulate(true),iterations_used_diagnostic(0)
{
    Print_Diagnostics();
    Print_Residuals(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLID_BODY_COLLECTION<TV>::
~SOLID_BODY_COLLECTION()
{
    delete &deformable_body_collection;
    delete &rigid_body_collection;
    delete &collision_body_list;
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Update_Simulated_Particles()
{
    rigid_body_collection.Update_Simulated_Particles();
    deformable_body_collection.Update_Simulated_Particles(*example_forces_and_velocities);
    int particles_number=deformable_body_collection.particles.Size();
    int rigid_particles_number=rigid_body_collection.rigid_body_particles.Size();

    ARRAY<bool> particle_is_simulated(particles_number);
    particle_is_simulated.Subset(deformable_body_collection.simulated_particles).Fill(true);

    ARRAY<bool> rigid_particle_is_simulated(rigid_particles_number);
    rigid_particle_is_simulated.Subset(rigid_body_collection.simulated_rigid_body_particles).Fill(true);
    for(int i=0;i<solid_force_collection.solids_forces.m;i++)
        solid_force_collection.solids_forces(i)->Update_Mpi(particle_is_simulated,rigid_particle_is_simulated,deformable_body_collection.mpi_solids);
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Update_Time_Varying_Material_Properties(const T time)
{
    example_forces_and_velocities->Update_Time_Varying_Material_Properties(time);
}
//#####################################################################
// Function Compute_Linear_Momentum 
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Compute_Linear_Momentum(TV& linear_momentum) const
{
    linear_momentum=TV();
    for(int i=0;i<deformable_body_collection.dynamic_particles.m;i++){int p=deformable_body_collection.dynamic_particles(i);
        linear_momentum+=deformable_body_collection.particles.mass(p)*deformable_body_collection.particles.V(p);}
    for(int i=0;i<rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=rigid_body_collection.dynamic_rigid_body_particles(i);
        linear_momentum+=rigid_body_collection.Rigid_Body(p).Mass()*rigid_body_collection.Rigid_Body(p).Twist().linear;}
}
//#####################################################################
// Function Compute_Energy
//#####################################################################
template<class TV> TV SOLID_BODY_COLLECTION<TV>::
Compute_Momentum() const
{
    TV momentum=deformable_body_collection.particles.V.Subset(deformable_body_collection.dynamic_particles).
        Weighted_Sum(deformable_body_collection.particles.mass.Subset(deformable_body_collection.dynamic_particles));
    for(int i=0;i<rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=rigid_body_collection.dynamic_rigid_body_particles(i);
        momentum+=rigid_body_collection.Rigid_Body(p).Momentum();}
    return momentum;
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Read(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const int static_frame,const bool include_static_variables,const bool read_rigid_body,const bool read_deformable_body,const bool read_from_every_process,ARRAY<int>* needs_init,ARRAY<int>* needs_destroy)
{
    std::string local_prefix=prefix;
    if(deformable_body_collection.mpi_solids && deformable_body_collection.mpi_solids->rank && !read_from_every_process){ // modify prefix to always read from the root's output
        std::string old_suffix="/"+STRING_UTILITIES::Value_To_String(deformable_body_collection.mpi_solids->rank+1),new_suffix="/1";
        size_t position=prefix.size()-old_suffix.size();
        PHYSBAM_ASSERT(prefix.substr(position)==old_suffix);
        local_prefix.replace(position,std::string::npos,new_suffix);}
    if(read_deformable_body){
        deformable_body_collection.Read(stream_type,local_prefix,local_prefix,frame,static_frame,include_static_variables,read_from_every_process);}
    if(read_rigid_body){
        rigid_body_collection.Read(stream_type,local_prefix,frame,needs_init,needs_destroy);}
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Write(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const int first_frame,const bool include_static_variables,const bool write_rigid_body,const bool write_deformable_body,const bool write_from_every_process,const bool output_interaction_pairs) const
{
    if(write_deformable_body){
        int static_frame=include_static_variables?frame:-1;
        bool write_static_variables=include_static_variables || frame==first_frame;
        deformable_body_collection.Write(stream_type,prefix,prefix,frame,static_frame,write_static_variables,write_from_every_process);
        ARRAY<FORCE_DATA<TV> > spring_data_list;
        for(int i=0;i<solid_force_collection.solids_forces.m;i++) solid_force_collection.solids_forces(i)->Add_Force_Data(spring_data_list);
        for(int i=0;i<deformable_body_collection.deformable_force_collection.deformables_forces.m;i++)
            deformable_body_collection.deformable_force_collection.deformables_forces(i)->Add_Force_Data(spring_data_list);
        std::string f=FILE_UTILITIES::Number_To_String(frame);
        if(spring_data_list.m!=0) FILE_UTILITIES::Write_To_File(stream_type,prefix+"/"+f+"/force_data",spring_data_list);}
    if(write_rigid_body)
        rigid_body_collection.Write(stream_type,prefix,frame);
    if(output_interaction_pairs)
        deformable_body_collection.triangle_repulsions.Output_Interaction_Pairs(stream_type,prefix+"/"+FILE_UTILITIES::Number_To_String(frame)+"/interaction_pairs");
}
//#####################################################################
namespace PhysBAM{
template class SOLID_BODY_COLLECTION<VECTOR<float,1> >;
template class SOLID_BODY_COLLECTION<VECTOR<float,2> >;
template class SOLID_BODY_COLLECTION<VECTOR<float,3> >;
template class SOLID_BODY_COLLECTION<VECTOR<double,1> >;
template class SOLID_BODY_COLLECTION<VECTOR<double,2> >;
template class SOLID_BODY_COLLECTION<VECTOR<double,3> >;
}

