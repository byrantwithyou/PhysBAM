//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// 1 one mesh, separate fragments
// 2 one mesh, unite fragments
// 3 multiple meshes, separate fragments
// 4 multiple meshes, unite fragments
// 4 multiple meshes, unite into 10 fragments
// 4 one mesh, unite into 10 fragments
//#####################################################################
#ifndef __FRAGMENT_TESTS__
#define __FRAGMENT_TESTS__
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>

namespace PhysBAM{

template<class T_input>
class FRAGMENT_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::stream_type;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes; // silence -Woverloaded-virtual

    GRID<TV> hair_layout_grid;
    bool unite_fragments;
    bool unite_some_fragments;
    bool use_single_mesh;


    FRAGMENT_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),unite_fragments(false),unite_some_fragments(false),use_single_mesh(false)
    {
    }

    // Unused callbacks

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    // deformable bodies
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    //RIGID_BODY_PARTICLES<TV>& rigid_body_particles=deformable_body_collection.rigid_body_particles;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    
    const int strands=25,strand_segments=50;
    //const int strands=2,strand_segments=50;
    hair_layout_grid.Initialize(TV_INT(strands,strands,strand_segments),RANGE<TV>(TV(-.5,0,0),TV(.5,1,2)));
    particles.Add_Elements(hair_layout_grid.Numbers_Of_Nodes().Product());
    
    int count=1;

    if(use_single_mesh){
        SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);deformable_body_collection.Add_Structure(&segmented_curve);
        for(int i=0;i<hair_layout_grid.counts.x;i++) for(int j=0;j<hair_layout_grid.counts.y;j++){
            for(int ij=1;ij<hair_layout_grid.counts.z;ij++){
                particles.X(count)=hair_layout_grid.X(TV_INT(i,j,ij));
                segmented_curve.mesh.elements.Append(VECTOR<int,2>(count,count+1));count++;}
            particles.X(count)=hair_layout_grid.X(TV_INT(i,j,hair_layout_grid.counts.z));count++;}
        T density=TV::dimension==1?1:TV::dimension==2?100:1000;
        SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(segmented_curve,density);}
    else{
        for(int i=0;i<hair_layout_grid.counts.x;i++) for(int j=0;j<hair_layout_grid.counts.y;j++){
            SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);deformable_body_collection.Add_Structure(&segmented_curve);
            for(int ij=1;ij<hair_layout_grid.counts.z;ij++){
                particles.X(count)=hair_layout_grid.X(TV_INT(i,j,ij));
                segmented_curve.mesh.elements.Append(VECTOR<int,2>(count,count+1));count++;}
            T density=TV::dimension==1?1:TV::dimension==2?100:1000;
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(segmented_curve,density);
            particles.X(count)=hair_layout_grid.X(TV_INT(i,j,hair_layout_grid.counts.z));count++;}}
    //RIGID_BODY<TV>& sphere_body=tests.Add_Rigid_Body("sphere",(T)1,(T)0.15);sphere_body.Frame().t=TV(0,(T)-1.2,0);

    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options()
{
    BASE::Register_Options();
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options()
{
    BASE::Parse_Options();
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    //SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    frame_rate=24;
    last_frame=(int)(64*frame_rate);
    solids_parameters.verbose_dt=true;
    solids_parameters.cfl=(T)4;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    if(test_number==4 || test_number==2) unite_fragments=true;
    if(test_number==1 || test_number==2 || test_number==6) use_single_mesh=true;
    if(test_number==5 || test_number==6) unite_some_fragments=true;
    last_frame=20;
    
    // geometry
    Get_Initial_Data();

    // make forces
    for(int i=0;i<deformable_body_collection.structures.m;i++){
        if(SEGMENTED_CURVE<TV>* curve=dynamic_cast<SEGMENTED_CURVE<TV>*>(deformable_body_collection.structures(i))){
            ARRAY<int>& referenced_nodes=*new ARRAY<int>; // hey craig, look a memory leak.  hi, andy, fix this one, too.
            curve->mesh.elements.Flattened().Get_Unique(referenced_nodes);
            solid_body_collection.solid_force_collection.Add_Force(Create_Edge_Springs(*curve,100/(1+sqrt((T)2)),(T)3));
            solid_body_collection.solid_force_collection.Add_Force(Create_Segment_Bending_Springs(*curve,100/(1+sqrt((T)2)),(T)3));
            solid_body_collection.solid_force_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,&referenced_nodes,0));}}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    for(int i=0;i<hair_layout_grid.counts.x;i++) for(int j=0;j<hair_layout_grid.counts.y;j++){
        int index=hair_layout_grid.counts.z*(j-1+(i-1)*hair_layout_grid.counts.y)+1;
        V(index)=TV();}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    for(int i=0;i<hair_layout_grid.counts.x;i++) for(int j=0;j<hair_layout_grid.counts.y;j++){
        int index=hair_layout_grid.counts.z*(j-1+(i-1)*hair_layout_grid.counts.y)+1;
        V(index)=TV();}
}
//#####################################################################
};
}
#endif
