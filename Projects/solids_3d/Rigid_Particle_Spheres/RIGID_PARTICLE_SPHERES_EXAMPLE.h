//#####################################################################
// Copyright 2006-2008, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_PARTICLE_SPHERES_EXAMPLE
//##################################################################### 
#ifndef __RIGID_PARTICLE_SPHERES_EXAMPLE__
#define __RIGID_PARTICLE_SPHERES_EXAMPLE__

#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/PARTICLE_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/IMPLICIT_ZERO_LENGTH_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class RIGID_PARTICLE_SPHERES_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
public:
    typedef T_input T;typedef VECTOR<T,3> TV;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::data_directory;
    using BASE::restart_frame;using BASE::output_directory;using BASE::solids_parameters;using BASE::write_last_frame;
    using BASE::fluids_parameters;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual

    SOLIDS_STANDARD_TESTS<TV> tests;
    SEGMENT_MESH segment_mesh;

    RIGID_PARTICLE_SPHERES_EXAMPLE(const STREAM_TYPE stream_type):
        BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection)
    {
    }

    // unused callbacks
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}

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
// Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    last_frame=1000;
    frame_rate=24;
    output_directory="Rigid_Particle_Spheres/output";
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-6;
    solids_parameters.implicit_solve_parameters.cg_iterations=500;
    segment_mesh.elements.Resize(1);
    solids_parameters.cfl=(T).5;

    switch(test_number){
      case 1:Deformable_Segment();break;
      case 2:Rigid_Particle_Segment();break;}

    deformable_body_collection.Add_Structure(new SEGMENTED_CURVE<TV>(segment_mesh,particles));

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

    // correct mass
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);


    LOG::cout<<"checking rigid body particle MASS"<<std::endl;
    for(int i=0;i<solid_body_collection.rigid_body_collection.rigid_body_particles.Size();i++){
        LOG::cout<<"solid_body_collection.rigid_body_collection.rigid_body_particles.mass("<<i<<")="<<solid_body_collection.rigid_body_collection.rigid_body_particles.mass(i)<<std::endl;
    }

    LOG::cout<<"checking bound particle effective  MASS"<<std::endl;
    for(int b=0;b<solid_body_collection.deformable_body_collection.binding_list.bindings.m;b++){
        int particle_index=solid_body_collection.deformable_body_collection.binding_list.bindings(b)->particle_index;
        LOG::cout<<"particles.mass.effective_mass("<<particle_index<<")="<<particles.effective_mass(particle_index)<<std::endl;
    }

    // add forces
    T overdamping_fraction=1;
    T stiffness=(T)1e4;
    IMPLICIT_ZERO_LENGTH_SPRINGS<TV>& implicit_zero_rest_length_springs=*Create_Edge_Zero_Length_Springs(deformable_body_collection.particles,segment_mesh,stiffness,overdamping_fraction);
/*    LINEAR_SPRINGS<TV>& implicit_zero_rest_length_springs=*Create_Edge_Springs(segment_mesh,deformable_body_collection.particles,stiffness,overdamping_fraction);
    ARRAY<T>::copy(0,implicit_zero_rest_length_springs.visual_restlength);implicit_zero_rest_length_springs.Clamp_Restlength(1);
    implicit_zero_rest_length_springs.Set_Overdamping_Fraction(overdamping_fraction);
    implicit_zero_rest_length_springs.use_implicit_velocity_independent_forces=true;*/

    solid_body_collection.solid_force_collection.Add_Force(&implicit_zero_rest_length_springs);
}
//#####################################################################
// Function Deformable_Segment
//#####################################################################
void Rigid_Particle_Segment()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    // add two rigid spheres
    for(int i=0;i<2;i++){
        RIGID_BODY<TV>& sphere=tests.Add_Rigid_Body("sphere",1,0);
        T sign=(T)(2*i-3);sphere.Frame().t+=TV::Axis_Vector(1)*(sign*2);
        segment_mesh.elements(1)(i)=particles.Add_Element();
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,segment_mesh.elements(1)(i),solid_body_collection.rigid_body_collection,sphere.particle_index,sign*TV()));}
}
//#####################################################################
// Function Deformable_Segment
//#####################################################################
void Deformable_Segment()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    for(int i=0;i<2;i++){
        T sign=(T)(2*i-3);
        int particle=particles.Add_Element();
        particles.X(particle)=TV::Axis_Vector(1)*(sign*2);
        particles.mass(particle)=1;
        int bound_particle=particles.Add_Element();
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new PARTICLE_BINDING<TV>(particles,bound_particle,particle));
        segment_mesh.elements(1)(i)=bound_particle;}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
}
};
}
#endif
