//#####################################################################
// Copyright 2004, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUASISTATIC_ATTACHMENT_EXAMPLE
//#####################################################################
#ifndef __QUASISTATIC_ATTACHMENT_EXAMPLE__
#define __QUASISTATIC_ATTACHMENT_EXAMPLE__

#include <Solids_And_Fluids/SOLIDS_3D_EXAMPLE.h>
namespace PhysBAM{

template <class T,class RW>
class QUASISTATIC_ATTACHMENT_EXAMPLE:public SOLIDS_3D_EXAMPLE<T,RW>
{
public:
    ARRAY<ARRAY<VECTOR_3D<T> > > target_positions,previous_positions,target_velocities,previous_velocities;
    ARRAY<SEGMENT_MESH> connections_to_quasistatic_particles;
    int current_frame;
    std::string quasistatic_output_directory;

    QUASISTATIC_ATTACHMENT_EXAMPLE()
    {
        last_frame=10*24;
        restart=false;restart_frame=0;   
        current_frame=restart_frame;
        cfl=(T)10;
        cg_tolerance=(T)1e-3;
        output_directory="Quasistatic_Attachment/output";
        quasistatic_output_directory="../quasistatics/Mattress/output";
        check_initial_mesh_for_self_intersection=false;
        perform_self_collision=false;
        solids_parameters.collide_with_interior=false;
    }

    ~QUASISTATIC_ATTACHMENT_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Deformable_Objects
//#####################################################################
void Initialize_Deformable_Objects()
{
    SOLIDS_3D_EXAMPLE<T,RW>::Initialize_Deformable_Objects();
    Add_Edge_Springs_To_Quasistatic_Particles(3000,(T).9);
}
//#####################################################################
// Function Add_Edge_Springs_To_Quasistatic_Particles
//#####################################################################
void Add_Edge_Springs_To_Quasistatic_Particles(const T stiffness=2e3,const T overdamping_fraction=1,const bool limit_time_step_by_strain_rate=true,
        const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,const T restlength_enlargement_fraction=0,const bool verbose=true)
{
    for(int i=0;i<solids_parameters.deformable_body_parameters.list.deformable_objects.m;i++){
    int number_quasistatic_particles=solids_parameters.deformable_body_parameters.list(i).tetrahedralized_volume->particles.array_collection->Size();
    ARRAY<int> links_to_quasistatic_particles(2,number_quasistatic_particles);
    for(int p=0;p<number_quasistatic_particles;p++){
        links_to_quasistatic_particles(1,p)=p;//quasistatic particle
        links_to_quasistatic_particles(2,p)=p+number_quasistatic_particles;}//dynamic particle
    connections_to_quasistatic_particles(1).Initizlize_Segment_Mesh(links_to_quasistatic_particles);//add link
    solids_parameters.deformable_body_parameters.list(i).Add_Edge_Springs(connections_to_quasistatic_particles(1),stiffness,overdamping_fraction,limit_time_step_by_strain_rate,
        max_strain_per_time_step,use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose);}
}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    //read quasistatic particles
    std::string prefix=quasistatic_output_directory+"/"+output_prefix;
    solids_parameters.deformable_body_parameters.list.template Read_Static_Variables<RW>(prefix);
    solids_parameters.deformable_body_parameters.list.template Read_Dynamic_Variables<RW>(prefix,restart_frame);

    //append new dynamic particles
    for(int i=0;i<solids_parameters.deformable_body_parameters.list.deformable_objects.m;i++){
        PARTICLES<T,VECTOR_3D<T> >& particles=solids_parameters.deformable_body_parameters.list(i).tetrahedralized_volume->particles;particles.Update_Velocity();
        particles.Increase_Array_Size(2*particles.array_collection->Size());
        for(int p=1;p<particles.array_collection->Size();p++)
        {int index=particles.array_collection->Add_Element();particles.X(index)=particles.X(p);particles.V(index)=particles.V(p);particles.mass(index)=particles.mass(p);}}
}
//#####################################################################
// Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time){
    int number_quasistatic_particles=solids_parameters.deformable_body_parameters.list(id_number).tetrahedralized_volume->particles.array_collection->Size()/2;
    for(int p=0;p<number_quasistatic_particles;p++) solids_parameters.deformable_body_parameters.list(id_number).tetrahedralized_volume->particles.V(p)=VECTOR_3D<T>(0,0,0);
}
//#####################################################################
// Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time){
    int number_quasistatic_particles=solids_parameters.deformable_body_parameters.list(id_number).tetrahedralized_volume->particles.array_collection->Size()/2;
    for(int p=0;p<number_quasistatic_particles;p++) solids_parameters.deformable_body_parameters.list(id_number).tetrahedralized_volume->particles.V(p)=VECTOR_3D<T>(0,0,0);
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time){
    T next_frame_time=initial_time+(current_frame+1-first_frame)/frame_rate;
    if(time>next_frame_time){
        current_frame++;
        for(int i=0;i<solids_parameters.deformable_body_parameters.list.deformable_objects.m;i++){
            for(int p=1;p<=solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume->particles.array_collection->Size()/2;p++)
        previous_positions(i)(p)=solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume->particles.X(p);
        previous_velocities(i)(p)=solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume->particles.V(p);}
    }
}
//#####################################################################
};
}
#endif
