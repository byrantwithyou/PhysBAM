//#####################################################################
// Copyright 2004, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "TP_EXAMPLE.h"
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NR0.h>

using namespace PhysBAM;

template<class T>
TP_EXAMPLE<T>::TP_EXAMPLE(int parameter) : ARTICULATED_RIGID_BODIES_3D_EXAMPLE<T,T>()
{
    last_frame=500;
    frame_rate=24*2;
    output_directory="TP/output";
    std::cout << "Frame rate: "<<frame_rate<<std::endl;
    arb.Set_Iterative_Tolerance((T)1e-6);
    //arb.Set_Extra_Iterations_Per_Contact_Level_Factor(100);
    //arb.Set_Extra_Iterations_Per_Shock_Propagation_Level_Factor(100);
    //arb.Set_Poststabilization_Iterations(100);
    write_last_frame=true;
}

template<class T>
TP_EXAMPLE<T>::~TP_EXAMPLE()
{}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void TP_EXAMPLE<T>::
Initialize_Rigid_Bodies()
{
    RANDOM_NR0 rg;
    int num_joints=0,num_bodies=0;
    int branch_joint=0,branch_body=0;
    VECTOR_3D<T> start_point,branch1,branch2,branch3,branch4;
    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(0,1,0),4,false,true,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(0,1,1),5,true,true,0,0,0,num_joints,num_bodies);
    branch1=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(1,-1,1),3,true,true,0,0,1,num_joints,num_bodies);branch_joint=num_joints-1;branch_body=num_bodies;
    start_point=Make_Lathe_Chain_2(branch1,VECTOR_3D<T>(0,-1,-1),2,true,true,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(1,1,-1),4,true,false,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(branch1,VECTOR_3D<T>(1,-1,1),3,true,true,branch_joint,branch_body,0,num_joints,num_bodies);

    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(0,-1,-1),5,true,true,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(-1,1,-1),4,true,true,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(-1,-.5,.5),6,true,true,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(0,0,1),5,true,true,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(1,1,0),6,true,true,0,0,0,num_joints,num_bodies);
    branch2=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(.5,-.5,-1),2,true,true,0,0,1,num_joints,num_bodies);branch_joint=num_joints-1;branch_body=num_bodies;
    start_point=Make_Lathe_Chain_2(branch2,VECTOR_3D<T>(1,-1,1),3,true,true,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(1,0,0),3,true,false,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(branch2,VECTOR_3D<T>(.5,-.5,-1),4,true,true,branch_joint,branch_body,0,num_joints,num_bodies);

    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(0,1,-.5),3,true,true,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(-1,0,0),6,true,true,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(0,0,1),3,true,true,0,0,0,num_joints,num_bodies);
    branch3=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(1,-.5,0),8,true,true,0,0,1,num_joints,num_bodies);branch_joint=num_joints-1;branch_body=num_bodies;
    start_point=Make_Lathe_Chain_2(branch3,VECTOR_3D<T>(.5,1,-1),6,true,false,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(branch3,VECTOR_3D<T>(1,-.5,0),2,true,true,branch_joint,branch_body,0,num_joints,num_bodies);

    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(.5,.5,-1),4,true,true,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(-1,.25,-.5),9,true,true,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(-1,-.5,1),6,true,true,0,0,0,num_joints,num_bodies);
    branch4=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(1,-1,.5),1,true,true,0,0,1,num_joints,num_bodies);branch_joint=num_joints-1;branch_body=num_bodies;
    start_point=Make_Lathe_Chain_2(branch4,VECTOR_3D<T>(0,.5,1),8,true,false,0,0,0,num_joints,num_bodies);
    start_point=Make_Lathe_Chain_2(branch4,VECTOR_3D<T>(1,-1,.5),4,true,true,branch_joint,branch_body,0,num_joints,num_bodies);

    start_point=Make_Lathe_Chain_2(start_point,VECTOR_3D<T>(1,1,0),4,true,false,0,0,0,num_joints,num_bodies);

    RIGID_BODY<TV>* rigid_body=0;
    rigid_body=Initialize_Rigid_Body("ground");
    rigid_body->position=VECTOR_3D<T>(0,-10,0);
    rigid_body->velocity=VECTOR_3D<T>(0,0,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
}
template<class T> VECTOR_3D<T> TP_EXAMPLE<T>::
Make_Lathe_Chain_2(VECTOR_3D<T> start_point,VECTOR_3D<T> direction,int number_of_links,bool start_joint,bool end_joint,int branch_joint,int branch_body,int branch_num,int& num_joints,int& num_bodies)
{
    direction.Normalize();
    RIGID_BODY<TV>* rigid_body=0;JOINT<TV>* joint;
    T link_length=4;
    VECTOR_3D<T> current_position=start_point;
    QUATERNION<T> orientation=QUATERNION<T>::Rotation_Quaternion(VECTOR_3D<T>(0,0,1),direction);// rotate the up vector to the direction vector
    for(int i=1;i<=number_of_links;i++){
        VECTOR_3D<T> link_position=current_position+direction*(link_length/2);
        current_position+=direction*link_length;
    
        if(i<number_of_links || end_joint){
            arb.joint_mesh.joint_description_list.Add_Joint_Description();
            joint=arb.joint_mesh.joints(arb.joint_mesh.Add_Joint());
            joint->joint_description.Set_Joint_To_Child_Frame(FRAME<T>(VECTOR_3D<T>(0,0,-2)));
            joint->joint_description.Set_Joint_To_Parent_Frame(FRAME<T>(VECTOR_3D<T>(0,0,2)));
        }

        rigid_body=Initialize_Rigid_Body("ARB/lathe_object");
        //rigid_body->position=start+VECTOR_3D<T>(2*sin(orient.x)*sin(orient.y),2*cos(orient.x),2*sin(orient.x)*cos(orient.y));
        rigid_body->position=link_position;
        rigid_body->orientation=orientation;
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(0.5);
    
    
    }
    if(start_joint) if(branch_joint) arb.Add_Articulation(branch_body,num_bodies+1,branch_joint);else arb.Add_Articulation(num_bodies,num_bodies+1,num_joints);
    for(int i=1;i<=number_of_links-1;i++) arb.Add_Articulation(num_bodies+i,num_bodies+i+1,num_joints+i);

    for(int j=1;j<=branch_num;j++){
        arb.joint_mesh.joint_description_list.Add_Joint_Description();
        joint=arb.joint_mesh.joints(arb.joint_mesh.Add_Joint());
        joint->joint_description.Set_Joint_To_Child_Frame(FRAME<T>(VECTOR_3D<T>(0,0,-2)));
        joint->joint_description.Set_Joint_To_Parent_Frame(FRAME<T>(VECTOR_3D<T>(0,0,2)));
    }

    num_joints+=number_of_links-1;
    if(end_joint) num_joints+=1;
    num_joints+=branch_num;
    num_bodies+=number_of_links;
    return current_position;
        
}
//#####################################################################
// Function Make_Lathe_Chain
//#####################################################################
template<class T> void TP_EXAMPLE<T>::
Make_Lathe_Chain(VECTOR_3D<T>& start,VECTOR_3D<T> orient,int num_in_chain,bool start_joint,bool end_joint,int& num_joints,int& num_bodies)
{
    RIGID_BODY<TV>* rigid_body=0;
    int current_body=1;
    JOINT<TV>* joint;
    while(current_body <= num_in_chain){

        if(current_body<num_in_chain || end_joint){
            arb.joint_mesh.joint_description_list.Add_Joint_Description();
            joint=arb.joint_mesh.joints(arb.joint_mesh.Add_Joint());
            joint->joint_description.Set_Joint_To_Child_Frame(FRAME<T>(VECTOR_3D<T>(0,0,-2)));
            joint->joint_description.Set_Joint_To_Parent_Frame(FRAME<T>(VECTOR_3D<T>(0,0,2)));
        }

        rigid_body=Initialize_Rigid_Body("ARB/lathe_object");
        //rigid_body->position=start+VECTOR_3D<T>(2*sin(orient.x)*sin(orient.y),2*cos(orient.x),2*sin(orient.x)*cos(orient.y));
        rigid_body->position=start+VECTOR_3D<T>(2*sin(orient.y)*cos(orient.x),2*sin(orient.x),2*cos(orient.x)*cos(orient.y));
        rigid_body->orientation=QUATERNION<T>(-orient.x,orient.y,orient.z);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(0.5);

        //start+=VECTOR_3D<T>(4*sin(orient.y)*sin(orient.x),4*cos(orient.x),4*sin(orient.x)*cos(orient.y));
        start+=VECTOR_3D<T>(4*sin(orient.y)*cos(orient.x),4*sin(orient.x),4*cos(orient.x)*cos(orient.y));
        current_body++;
    }

    if(start_joint) arb.Add_Articulation(num_bodies,num_bodies+1,num_joints);
    for(int i=1;i<=num_in_chain-1;i++) arb.Add_Articulation(num_bodies+i,num_bodies+i+1,num_joints+i);

    num_joints+=num_in_chain-1;
    if(end_joint) num_joints+=1;
    num_bodies+=num_in_chain;
}

namespace PhysBAM
{
    template TP_EXAMPLE<double>;
    template TP_EXAMPLE<float>;
}
