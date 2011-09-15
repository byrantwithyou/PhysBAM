//#####################################################################
// Copyright 2004-2008, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_PARTICLE_MESH_EXAMPLE
//##################################################################### 
#ifndef __RIGID_PARTICLE_MESH_EXAMPLE__
#define __RIGID_PARTICLE_MESH_EXAMPLE__

#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class RIGID_PARTICLE_MESH_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
public:
    typedef T_input T;
    typedef VECTOR<T,3> TV;typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::data_directory;
    using BASE::restart_frame;using BASE::output_directory;using BASE::solids_parameters;using BASE::write_last_frame;
    using BASE::fluids_parameters;using BASE::solid_body_collection;using BASE::parse_args;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes; // silence -Woverloaded-virtual

    SOLIDS_STANDARD_TESTS<TV> tests;
    SEGMENT_MESH segment_mesh;
    int num_rows,num_cols;
    T stiffness,overdamping_fraction;

    RIGID_PARTICLE_MESH_EXAMPLE(const STREAM_TYPE stream_type):
        BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),num_rows(8),num_cols(8),stiffness(1e4),overdamping_fraction(1)
    {
    }

    // unused callbacks
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
void Add_Rigid_Body(const std::string& rigid_body_name,const std::string& filename,const FRAME<TV>& frame,const bool with_phi)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body(filename,(T).15,0,with_phi);
    rigid_body.X()=frame.t;
    rigid_body.Rotation()=frame.r;
    rigid_body.Set_Coefficient_Of_Restitution(0);
    rigid_body.Set_Name(rigid_body_name);
    int gravity_particle=particles.array_collection->Add_Element();
    particles.mass(gravity_particle)=rigid_body.Mass();
    solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,gravity_particle,solid_body_collection.rigid_body_collection,rigid_body.particle_index,TV()));
}
//#####################################################################
// Function Add_Point_Joint
//#####################################################################
void Add_Point_Joint(const int joint_id,const int plank_id,const TV& plank_object_space_position,const int sphere_id,const TV& sphere_object_space_position)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    // add deformable particles for each of the sphere and plank
    int sphere_particle=particles.array_collection->Add_Element();
    int plank_particle=particles.array_collection->Add_Element();
    segment_mesh.elements.Append(VECTOR<int,2>(sphere_particle,plank_particle));
    
    // add bindings to their rigid body particles
    solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,sphere_particle,solid_body_collection.rigid_body_collection,sphere_id,sphere_object_space_position));
    solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,plank_particle,solid_body_collection.rigid_body_collection,plank_id,plank_object_space_position));
}
//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options()
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-rows",1,"","num rows, num cols");
    parse_args->Add_Double_Argument("-stiffen",1,"","stiffness multiplier for various tests");
    parse_args->Add_Double_Argument("-dampen",1,"","damping multiplier for various tests");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options()
{
    BASE::Parse_Options();
    if(parse_args->Is_Value_Set("-rows")) num_rows=num_cols=parse_args->Get_Integer_Value("-rows");
    if(parse_args->Is_Value_Set("-dampen")) overdamping_fraction=(T)parse_args->Get_Double_Value("-dampen");
    if(parse_args->Is_Value_Set("-stiffen")) stiffness=(T)parse_args->Get_Double_Value("-stiffen");
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.cfl=(T).2;
    last_frame=10000;
    frame_rate=24;
    output_directory="Rigid_Particle_Mesh/output";
    LOG::cout<<"Frame rate: "<<frame_rate<<std::endl;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-6;
    solids_parameters.implicit_solve_parameters.cg_iterations=500;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;

    T x_shift=-8,y_shift=0,z_shift=-8;

    bool with_phi=false;

    for(int row=0;row<num_rows;row++){
        for(int col=0;col<num_cols;col++) Add_Rigid_Body("mesh","plank",FRAME<TV>(TV(x_shift+1+2*col,y_shift,z_shift+2*row),ROTATION<TV>((T)pi/2,TV(0,1,0))),with_phi);
        for(int col=0;col<num_cols+1;col++) Add_Rigid_Body("mesh_joint","sphere",FRAME<TV>(TV(x_shift+2*col,y_shift,z_shift+2*row)),with_phi);
        for(int col=0;col<num_cols+1;col++) Add_Rigid_Body("mesh","plank",FRAME<TV>(TV(x_shift+2*col,y_shift,z_shift+1+2*row)),with_phi);
        if(row==num_rows-1) for(int col=0;col<num_cols+1;col++) Add_Rigid_Body("mesh_joint","sphere",FRAME<TV>(TV(x_shift+2*col,y_shift,z_shift+2+2*row)),with_phi);}
    for(int col=0;col<num_cols;col++) Add_Rigid_Body("mesh","plank",FRAME<TV>(TV(x_shift+1+2*col,y_shift,z_shift+2*num_rows),ROTATION<TV>((T)pi/2,TV(0,1,0))),with_phi);

    int i=0;
    for(int row=0;row<num_rows+1;row++) for(int col=0;col<num_cols+1;col++){
        if(row==0){
            if(col==0){
                Add_Point_Joint(++i,1,TV(0,0,-1),num_cols+1,TV());
                Add_Point_Joint(++i,2*num_cols+2,TV(0,0,-1),num_cols+1,TV());}
            else if(col==num_cols){
                Add_Point_Joint(++i,1,TV(0,0,-1),num_cols+1,TV());
                Add_Point_Joint(++i,num_cols,TV(0,0,1),2*num_cols+1,TV());
                Add_Point_Joint(++i,3*num_cols+2,TV(0,0,-1),2*num_cols+1,TV());}
            else{
                Add_Point_Joint(++i,col,TV(0,0,1),num_cols+1+col,TV());
                Add_Point_Joint(++i,col+1,TV(0,0,-1),num_cols+1+col,TV());
                Add_Point_Joint(++i,2*num_cols+2+col,TV(0,0,-1),num_cols+1+col,TV());}}
        else if(row==num_rows){
            int total=num_cols*num_rows+2*(num_cols+1)*num_rows;
            if(col==0){
                Add_Point_Joint(++i,total-num_cols,TV(0,0,1),total+1,TV());
                Add_Point_Joint(++i,total+num_cols+2,TV(0,0,-1),total+1,TV());}
            else if(col==num_cols){
                Add_Point_Joint(++i,total,TV(0,0,1),total+num_cols+1,TV());
                Add_Point_Joint(++i,total+2*num_cols+1,TV(0,0,1),total+num_cols+1,TV());}
            else{
                Add_Point_Joint(++i,total-num_cols+col,TV(0,0,1),total+1+col,TV());
                Add_Point_Joint(++i,total+num_cols+1+col,TV(0,0,1),total+1+col,TV());
                Add_Point_Joint(++i,total+num_cols+2+col,TV(0,0,-1),total+1+col,TV());}}
        else{
            int total=num_cols*row+2*(num_cols+1)*row;
            if(col==0){
                Add_Point_Joint(++i,total+1,TV(0,0,-1),total+num_cols+1,TV());
                Add_Point_Joint(++i,total-num_cols,TV(0,0,1),total+num_cols+1,TV());
                Add_Point_Joint(++i,total+2*num_cols+2,TV(0,0,-1),total+num_cols+1,TV());
            }
            else if(col==num_cols){
                Add_Point_Joint(++i,total,TV(0,0,1),total+2*num_cols+1,TV());
                Add_Point_Joint(++i,total+num_cols,TV(0,0,1),total+2*num_cols+1,TV());
                Add_Point_Joint(++i,total+3*num_cols+2,TV(0,0,-1),total+2*num_cols+1,TV());
            }
            else{
                Add_Point_Joint(++i,total+col,TV(0,0,1),total+num_cols+1+col,TV());
                Add_Point_Joint(++i,total+col+1,TV(0,0,-1),total+num_cols+1+col,TV());
                Add_Point_Joint(++i,total-num_cols+col,TV(0,0,1),total+num_cols+1+col,TV());
                Add_Point_Joint(++i,total+2*num_cols+2+col,TV(0,0,-1),total+num_cols+1+col,TV());}}}

    deformable_body_collection.deformable_geometry.Add_Structure(new SEGMENTED_CURVE<TV>(segment_mesh,particles));

#if 0
    int id=0;
    RIGID_BODY<TV>* rigid_body=0;

    // boxes
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    rigid_body=parameters.rigid_body_parameters.list.rigid_bodies(id);
    rigid_body->frame.t=TV(-3,(T)1.25,-3);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    rigid_body=parameters.rigid_body_parameters.list.rigid_bodies(id);
    rigid_body->frame.t=TV(-3,(T)1.25,3);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    rigid_body=parameters.rigid_body_parameters.list.rigid_bodies(id);
    rigid_body->frame.t=TV(3,(T)1.25,3);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box");
    rigid_body=parameters.rigid_body_parameters.list.rigid_bodies(id);
    rigid_body->frame.t=TV(3,(T)1.25,-3);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
#endif

    // correct number nodes
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);

    // add forces
    solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));

    IMPLICIT_ZERO_LENGTH_SPRINGS<TV>& zero_length_springs=*Create_Edge_Zero_Length_Springs(deformable_body_collection.particles,segment_mesh,stiffness,overdamping_fraction);
/*    LINEAR_SPRINGS<TV>& zero_length_springs=*Create_Edge_Springs(segment_mesh,deformable_body_collection.particles,stiffness,overdamping_fraction);
    ARRAY<T>::copy(0,zero_length_springs.visual_restlength);zero_length_springs.Clamp_Restlength(1);
    zero_length_springs.Set_Overdamping_Fraction(overdamping_fraction);
    zero_length_springs.use_implicit_velocity_independent_forces=true;*/

    solid_body_collection.Add_Force(&zero_length_springs);

//    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    solid_body_collection.Update_Simulated_Particles();
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    twist(num_cols+1)=TWIST<TV>();
    twist(2*num_cols+1)=TWIST<TV>();
    twist(num_cols*num_rows+2*(num_cols+1)*num_rows+1)=TWIST<TV>();
    twist(num_cols*num_rows+2*(num_cols+1)*num_rows+num_cols+1)=TWIST<TV>();
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    twist(num_cols+1)=TWIST<TV>();
    twist(2*num_cols+1)=TWIST<TV>();
    twist(num_cols*num_rows+2*(num_cols+1)*num_rows+1)=TWIST<TV>();
    twist(num_cols*num_rows+2*(num_cols+1)*num_rows+num_cols+1)=TWIST<TV>();
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{}
//#####################################################################
};
}
#endif
