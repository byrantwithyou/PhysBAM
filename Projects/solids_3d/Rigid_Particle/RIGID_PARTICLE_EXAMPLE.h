//#####################################################################
// Copyright 2006-2008, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_PARTICLE_EXAMPLE
//#####################################################################
#ifndef __RIGID_PARTICLE_EXAMPLE__
#define __RIGID_PARTICLE_EXAMPLE__

#include <Core/Log/DEBUG_UTILITIES.h>
#include <Deformables/Bindings/RIGID_BODY_BINDING.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Deformables/Forces/SEGMENT_BENDING_ELEMENTS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
namespace PhysBAM{

template<class T_input>
class RIGID_PARTICLE_EXAMPLE:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
public:
    typedef T_input T;typedef VECTOR<T,3> TV;
    typedef SOLIDS_EXAMPLE<TV> BASE;
    typedef typename TV::SPIN T_SPIN;
    using BASE::last_frame;using BASE::restart;using BASE::viewer_dir;using BASE::solids_parameters;using BASE::data_directory;
    using BASE::stream_type;using BASE::solid_body_collection;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes; // silence -Woverloaded-virtual
    using BASE::user_last_frame;
    
    T initial_height;
    bool use_ground;
    int ground_id;

    RIGID_PARTICLE_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args),initial_height((T)5),use_ground(true)
    {
        parse_args.Parse();
    }

    ~RIGID_PARTICLE_EXAMPLE()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
virtual void Get_Initial_Data()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);

    // add particles and set up segmented curve
    for(int i=0;i<20;i++){
        particles.Add_Element();
        particles.mass(i)=1;
        particles.X(i)=TV((T)i,initial_height,0);
        if(i>1) segmented_curve.mesh.elements.Append(VECTOR<int,2>(i,i-1));}
    deformable_body_collection.Add_Structure(&segmented_curve);

    // rigid bodies
    for(int i=0;i<3;i++){
        solid_body_collection.rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies/subdivided_box",(T)1.0);
        rigid_body_collection.rigid_body_particles.frame(i).t=TV((T)5*i+(T).5,initial_height,0);
        rigid_body_collection.rigid_body_particles.frame(i).r=ROTATION<TV>::From_Rotation_Vector(TV());
        rigid_body_collection.Rigid_Body(i).coefficient_of_friction=0;
        rigid_body_collection.Rigid_Body(i).Set_Coefficient_Of_Restitution((T).5);
        T mass_scale_factor=10/rigid_body_collection.rigid_body_particles.mass(i);
        rigid_body_collection.Rigid_Body(i).Set_Rigid_Mass(rigid_body_collection.Rigid_Body(i).Rigid_Mass()*mass_scale_factor);
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,5*i,rigid_body_collection,i,TV((T)-.5,0,0)));
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,5*i+1,rigid_body_collection,i,TV((T).5,0,0)));}

    segmented_curve.Update_Number_Nodes();
    if(use_ground){
        ground_id=solid_body_collection.rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies/ground",(T).045);
        solid_body_collection.rigid_body_collection.Rigid_Body(ground_id).coefficient_of_friction=(T).2;
        solid_body_collection.rigid_body_collection.Rigid_Body(ground_id).Set_Coefficient_Of_Restitution((T).5);
        solid_body_collection.rigid_body_collection.Rigid_Body(ground_id).Is_Kinematic()=true;}

    // correct mass
    deformable_body_collection.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);

    // collisions
//    solids_parameters.collision_body_list.Add_Bodies(solid_body_collection.rigid_body_collection);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    if(!user_last_frame) last_frame=2400;
    restart=0;  
    solids_parameters.cfl=(T)1;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
    if(!this->user_output_directory)
        viewer_dir.output_directory="Rigid_Particle/output";
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;

    Get_Initial_Data();

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    SEGMENTED_CURVE<TV>& segmented_curve=deformable_body_collection.template Find_Structure<SEGMENTED_CURVE<TV>&>();

    solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));

    solid_body_collection.Add_Force(Create_Edge_Springs(segmented_curve,(T)1e2,(T)1));
    //solid_body_collection.Add_Force(Create_Bending_Elements(segmented_curve,(T)1e1,(T)1));

    solid_body_collection.Update_Simulated_Particles();

    SOLIDS_EXAMPLE<TV>::Initialize_Bodies();
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override
{
    if(use_ground) twist(Value(ground_id))=TWIST<TV>();
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override
{
    if(use_ground) twist(Value(ground_id))=TWIST<TV>();
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) override
{
    if(id==ground_id) frame.t.x=solid_body_collection.rigid_body_collection.rigid_body_particles.frame(2).t.x;
}
//#####################################################################
};
}
#endif
