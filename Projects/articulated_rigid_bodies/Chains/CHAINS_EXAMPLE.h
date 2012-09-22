//#####################################################################
// Copyright 2004-2008, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CHAINS_EXAMPLE
//##################################################################### 
#ifndef __CHAINS_EXAMPLE__
#define __CHAINS_EXAMPLE__

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/RIGID_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class CHAINS_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::fluids_parameters;using BASE::output_directory;using BASE::solids_parameters;using BASE::write_last_frame;using BASE::data_directory;using BASE::stream_type;
    using BASE::restart;using BASE::initial_time;using BASE::first_frame;using BASE::last_frame;using BASE::restart_frame;using BASE::frame_rate;using BASE::solid_body_collection;
    using BASE::test_number;using BASE::parse_args;using BASE::Set_External_Velocities; // silence -Woverloaded-virtual
    
    ARTICULATED_RIGID_BODY<TV>* arb;
    SOLIDS_STANDARD_TESTS<TV> tests;
    RIGID_BODY<TV> *square1,*square2;
    int num_poles,num_rings,selection;
    bool use_rigid_deformable_evolution_old;

    CHAINS_EXAMPLE(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),arb(0),tests(*this,solid_body_collection),square1(0),square2(0),num_poles(5),num_rings(20),selection(0),
        use_rigid_deformable_evolution_old(false)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        //solids_parameters.rigid_body_parameters.particle_partition_size=6;
        solids_parameters.cfl=(T).1;
        //solids_parameters.perform_collision_body_collisions=false;


        last_frame=500;
        frame_rate=24*2;
        write_last_frame=true;

        //artificial_maximum_speed=30;
        std::cout<<"Frame rate: "<<frame_rate<<std::endl;
    }
    
    virtual ~CHAINS_EXAMPLE()
    {
        delete &arb;
    }

    // Unused callbacks
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}

    void Register_Options() PHYSBAM_OVERRIDE
    {
        BASE::Register_Options();
        parse_args->Add("-selection",&selection,"value","selection");
    }
    void Parse_Options() PHYSBAM_OVERRIDE
    {
        BASE::Parse_Options();
        output_directory="Chains/output";output_directory+=selection==0?"_blocks_chain":"_lathe_chains";
    }
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    int num_joints=0,num_bodies=0;
    T up_shift=60;
    RANDOM_NUMBERS<T> random_generator;

    arb=&solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb->Set_Iterative_Tolerance((T)1e-6);
    //arb->Set_Extra_Iterations_Per_Contact_Level_Factor(100);
    //arb->Set_Extra_Iterations_Per_Shock_Propagation_Level_Factor(100);
    //arb->Set_Poststabilization_Iterations(100);
    arb->Set_Do_Final_Pass(false);
    //arb->Set_Use_Shock_Propagation(false);
    std::cout<<"Use shock propagation: "<<arb->use_shock_propagation<<std::endl;

    if(selection==0) Make_Block_Chain(TV(),ROTATION<TV>(),num_joints,num_bodies);
    else{
        // rings
        for(int j=0;j<num_rings;j++) for(int i=-1;i<=1;i++) for(int k=-1;k<=1;k++){
            T rand=random_generator.Get_Number()*2-1;
            if(j%2) tests.Make_Lathe_Chain(FRAME<TV>(TV((T)i*10,up_shift+j*10,k*4+rand))/*,"chains1"*/); // random x
            else tests.Make_Lathe_Chain(FRAME<TV>(TV(i*4+rand,up_shift+j*10,(T)k*10),ROTATION<TV>((T)pi/2,TV(0,1,0)))/*,"chains2"*/);} // random y
        // poles
        for(int i=0;i<num_poles;i++) for(int j=0;j<num_poles;j++){
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("Rings_Test/medium_cylinder",1,(T).6);
            rigid_body.is_static=true;
            rigid_body.name=STRING_UTILITIES::string_sprintf("pole %d %d",i,j);
            rigid_body.Frame().t=TV((i-(num_poles+1)/(T)2)*7,10,(j-(num_poles+1)/(T)2)*7);}}

    tests.Add_Ground((T).5,0,(T).5);
    tests.Add_Gravity();

    solid_body_collection.Update_Simulated_Particles();

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Initialize_Bodies();
}
//#####################################################################
// Function Make_Block_Chain
//#####################################################################
void Make_Block_Chain(TV shift,ROTATION<TV> orient,int& num_joints,int& num_bodies)
{
    RIGID_BODY<TV>* rigid_bodies[13];
    JOINT<TV>* joints[12];
    TV positions[12]={TV(-2,4,0),TV(0,4,0),TV(2,4,0),TV(4,2,0),TV(4,0,0),TV(4,-2,0),TV(2,-4,0),TV(0,-4,0),TV(-2,-4,0),TV(-4,-2,0),TV(-4,0,0),TV(-4,2,0)};
    for(int i=0;i<12;i++){
        int id=solid_body_collection.rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/subdivided_box");
        rigid_bodies[i]=&arb->rigid_body_collection.Rigid_Body(id);
        rigid_bodies[i]->Frame().t=orient.Rotate(positions[i])+shift;
        rigid_bodies[i]->Frame().r=orient;
        rigid_bodies[i]->Set_Coefficient_Of_Restitution(0.5);}

    rigid_bodies[12]=rigid_bodies[0];
    for(int i=0;i<12;i++){joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(rigid_bodies[i]->particle_index,rigid_bodies[i+1]->particle_index,joints[i]);}

    for(int j=0;j<3;j++) joints[j]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,-1,1)));
    for(int j=3;j<6;j++) joints[j]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,-1,1)));
    for(int j=6;j<9;j++) joints[j]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,1)));
    for(int j=9;j<12;j++) joints[j]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));
    for(int j=0;j<2;j++) joints[j]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,1)));
    for(int j=2;j<5;j++) joints[j]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,1,1)));
    for(int j=5;j<8;j++) joints[j]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,1,1)));
    for(int j=8;j<12;j++) joints[j]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,1)));
    joints[11]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,1)));
    num_joints+=12;num_bodies+=12;

    //set velocities
    arb->rigid_body_collection.rigid_body_particle.angular_momentum(3)=TV(10,0,0);
    arb->rigid_body_collection.rigid_body_particle.angular_momentum(9)=TV(-10,20,0);
    arb->rigid_body_collection.rigid_body_particle.angular_momentum(5)=TV(10,10,0);
    arb->rigid_body_collection.rigid_body_particle.angular_momentum(2)=TV(5,0,10);
    arb->rigid_body_collection.rigid_body_particle.angular_momentum(1)=TV(0,0,0);
    arb->rigid_body_collection.rigid_body_particle.angular_momentum(8)=TV(2,0,10);
}
//#####################################################################
};
}
#endif
