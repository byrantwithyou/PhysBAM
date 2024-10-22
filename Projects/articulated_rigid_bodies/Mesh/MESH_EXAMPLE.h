//#####################################################################
// Copyright 2004-2008, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MESH_EXAMPLE
//##################################################################### 
#ifndef __MESH_EXAMPLE__
#define __MESH_EXAMPLE__

#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Joints/POINT_JOINT.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include "../ARB_PARAMETERS.h"
namespace PhysBAM{

template<class T_input>
class MESH_EXAMPLE:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
public:
    typedef T_input T;
    typedef VECTOR<T,3> TV;typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::viewer_dir;using BASE::solids_parameters;
    using BASE::data_directory;using BASE::stream_type;using BASE::solid_body_collection;
    using BASE::Set_External_Velocities; // silence -Woverloaded-virtual
    using BASE::user_last_frame;
    
    ARTICULATED_RIGID_BODY<TV>* arb;
    RIGID_BODY<TV> *handle[4];
    TV pos[4];
    T increment;
    PARAMETER_LIST parameter_list;
    SOLIDS_STANDARD_TESTS<TV> tests;

    MESH_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args):
        BASE(stream_type_input,parse_args),tests(stream_type_input,data_directory,solid_body_collection)
    {
        parse_args.Parse();
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=(T).1;

        if(!user_last_frame) last_frame=1000;
        if(!this->user_frame_rate) frame_rate=24;
        std::cout<<"Frame rate: "<<frame_rate<<std::endl;
        increment=0;
        tests.data_directory=data_directory;
        if(!this->user_output_directory)
            viewer_dir.output_directory="Mesh/output";
    }

    virtual ~MESH_EXAMPLE()
    {
        delete arb;
    }
//#####################################################################
// Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    int id(0);
    T x_shift=-8,y_shift=0,z_shift=-8;
    RIGID_BODY<TV>* rigid_body;
    int num_rows=8,num_cols=8;

    bool with_phi=false;

    arb=&solid_body_collection.rigid_body_collection.articulated_rigid_body;
    //arb->Set_Extra_Iterations_Per_Contact_Level_Factor(2000);
    //arb->Set_Extra_Iterations_Per_Shock_Propagation_Level_Factor(2000);
    //arb->Set_Poststabilization_Iterations(2000);
    //arb->Use_Epsilon_Scale(false);
    //solids_parameters.rigid_body_collision_parameters.use_epsilon_scaling=false;
    //solids_parameters.rigid_body_collision_parameters.use_epsilon_scaling_for_level=false;
    //arb->post_stabilization_in_pre_stabilization_configuration=true;
    //arb->epsilon_scaling_for_post_stabilization=true;
    arb->Set_Iterative_Tolerance((T)1e-5);
    arb->Set_Contact_Level_Iterations(10);
    arb->Set_Shock_Propagation_Level_Iterations(10);
    arb->Set_Poststabilization_Iterations(10);
    arb->Set_Use_Shock_Propagation(false);
    arb->Set_Do_Final_Pass(false);

    ARB_PARAMETERS::Read_Common_Parameters("Mesh/example.param",*this,parameter_list);

    for(int row=0;row<num_rows;row++){
        for(int col=0;col<num_cols;col++){
            id=rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies/"+"plank",(T).15,true,with_phi);
            rigid_body=&arb->rigid_body_collection.Rigid_Body(id);
            rigid_body->Frame().t=TV(x_shift+1+2*col,y_shift,z_shift+2*row);
            rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0));
            rigid_body->Set_Coefficient_Of_Restitution(0);
            rigid_body->name="mesh";}
        for(int col=0;col<num_cols+1;col++){
            id=rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies/"+"sphere",(T).15,true,with_phi);
            rigid_body=&arb->rigid_body_collection.Rigid_Body(id);
            rigid_body->Frame().t=TV(x_shift+2*col,y_shift,z_shift+2*row);
            rigid_body->Set_Coefficient_Of_Restitution(0);
            rigid_body->name="mesh_joint";}
        for(int col=0;col<num_cols+1;col++){
            id=rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies/"+"plank",(T).15,true,with_phi);
            rigid_body=&arb->rigid_body_collection.Rigid_Body(id);
            rigid_body->Frame().t=TV(x_shift+2*col,y_shift,z_shift+1+2*row);
            rigid_body->Set_Coefficient_Of_Restitution(0);
            rigid_body->name="mesh";}            
        if(row==num_rows-1){
            for(int col=0;col<num_cols+1;col++){
                id=rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies/"+"sphere",(T).15,true,with_phi);
                rigid_body=&arb->rigid_body_collection.Rigid_Body(id);
                rigid_body->Frame().t=TV(x_shift+2*col,y_shift,z_shift+2+2*row);
                rigid_body->Set_Coefficient_Of_Restitution(0);
                rigid_body->name="mesh_joint";}}}
    for(int col=0;col<num_cols;col++){
        id=rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies/"+"plank",(T).15,true,with_phi);
        rigid_body=&arb->rigid_body_collection.Rigid_Body(id);
        rigid_body->Frame().t=TV(x_shift+1+2*col,y_shift,z_shift+2*num_rows);
        rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0));
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->name="mesh";}

    int num_joints=2+3*(num_cols-1) + 2+3*(num_cols-1) + 6*(num_rows-1)+4*num_cols*(num_rows-1);
    JOINT<TV>** joints;joints=new JOINT<TV>*[num_joints];
    int i=0;

    for(int row=0;row<num_rows+1;row++){
        for(int col=0;col<num_cols+1;col++){
            if(row==0){
                if(col==0){
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_cols+1),int(1),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,-1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_cols+1),int(2*num_cols+2),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,-1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))));}
               else if(col==num_cols){
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(2*num_cols+1),int(num_cols),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(2*num_cols+1),int(3*num_cols+2),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,-1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))));}
                else{
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_cols+1+col),int(col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_cols+1+col),int(col+1),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,-1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_cols+1+col),int(2*num_cols+2+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,-1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))));}
            }
            else if(row==num_rows){
                int total=num_cols*num_rows+2*(num_cols+1)*num_rows;
                if(col==0){
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+1),int(total-num_cols),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))));

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+1),int(total+num_cols+2),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,-1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,1,0))));} // --
                else if(col==num_cols){
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1),int(total),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))));

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1),int(total+2*num_cols+1),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,1,0))));} // --
                else{
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+1+col),int(total-num_cols+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))));

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+1+col),int(total+num_cols+1+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,1,0)))); // --
 
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+1+col),int(total+num_cols+2+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,-1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,1,0))));}} // --
            else{
                int total=num_cols*row+2*(num_cols+1)*row;
                if(col==0){
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1),int(total+1),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,-1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1),int(total-num_cols),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))));

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1),int(total+2*num_cols+2),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,-1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))));}
                else if(col==num_cols){
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+2*num_cols+1),int(total),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))));

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+2*num_cols+1),int(total+num_cols),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+2*num_cols+1),int(total+3*num_cols+2),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,-1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))));}
                else{
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1+col),int(total+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1+col),int(total+col+1),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,-1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1+col),int(total-num_cols+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))));
 
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1+col),int(total+2*num_cols+2+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,-1),ROTATION<TV>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))));}}}}

    //handles
    handle[0]=&arb->rigid_body_collection.Rigid_Body(num_cols+1);
    handle[1]=&arb->rigid_body_collection.Rigid_Body(2*num_cols+1);
    handle[2]=&arb->rigid_body_collection.Rigid_Body(num_cols*num_rows+2*(num_cols+1)*num_rows+1);
    handle[3]=&arb->rigid_body_collection.Rigid_Body(num_cols*num_rows+2*(num_cols+1)*num_rows+num_cols+1);
    for(int i=0;i<4;i++){handle[i]->is_static=true;pos[i]=handle[i]->Frame().t;}

#if 1
    // boxes
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",1,.5);
    rigid_body->Frame().t=TV(-3,1.25,-3);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);

    rigid_body=&tests.Add_Rigid_Body("subdivided_box",1,.5);
    rigid_body->Frame().t=TV(-3,1.25,3);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);

    rigid_body=&tests.Add_Rigid_Body("subdivided_box",1,.5);
    rigid_body->Frame().t=TV(3,1.25,3);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);

    rigid_body=&tests.Add_Rigid_Body("subdivided_box",1,.5);
    rigid_body->Frame().t=TV(3,1.25,-3);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
#endif

    // the ground!
    tests.Add_Ground((T).3,-10,(T).5);
    tests.Add_Gravity();
    // TODO: this example also add ether drag, but when basic forces were removed we didn't add it.
}

//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) override
{
     
    TV center=(T).25*(pos[0]+pos[1]+pos[2]+pos[3]),top(center.x,15,center.z);
    if(frame<12){
        T increment=frame*(T).01;
        for(int i=0;i<4;i++) handle[i]->Frame().t=TV((1-increment)*pos[i]);}
    
}
};
}
#endif
