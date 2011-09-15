//#####################################################################
// Copyright 2004, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include "SNAKE_EXAMPLE.h"

using namespace PhysBAM;

template<class T,class RW>
SNAKE_EXAMPLE<T,RW>::SNAKE_EXAMPLE() : SOLIDS_FLUIDS_EXAMPLE_3D<RW>(FLUIDS_PARAMETERS_3D<T>::NONE)
{
    solids_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;

    last_frame=500;
    frame_rate=24;
    output_directory="Snake/output";
    solids_parameters.rigid_body_parameters.artificial_maximum_speed=20;
    std::cout << "Frame rate: "<<frame_rate<<std::endl;

    arb=new ARTICULATED_RIGID_BODY<TV>(this->solids_parameters.rigid_body_parameters.list);
    this->solids_parameters.rigid_body_parameters.Set_Articulated_Rigid_Body(arb);

    arb->Set_Iterative_Tolerance((T)1e-6);
    //arb->Set_Extra_Iterations_Per_Contact_Level_Factor(100);
    //arb->Set_Extra_Iterations_Per_Shock_Propagation_Level_Factor(100);
    //arb->Set_Poststabilization_Iterations(100);
    arb->Set_Do_Final_Pass(false);
    write_last_frame=true;
    //particle_partition_size=6;
    current_frame=0;
}

template<class T,class RW>
SNAKE_EXAMPLE<T,RW>::~SNAKE_EXAMPLE()
{
    delete arb;
}

//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T,class RW> void SNAKE_EXAMPLE<T,RW>::
Initialize_Bodies()
{
//    int num_joints=0,num_bodies=0;
    
    RIGID_BODY<TV> *rigid_body = 0;int id=0;
    // make a center standing pole to connect joint to
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"Rings_Test/cylinder_revolve");
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body->frame.t=VECTOR_3D<T>(0,0,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Name("pole");
    rigid_body->is_static=true;
    
//REPLICATE THIS LINE TO ADD MORE SNAKES, JUST GIVE THEM A DIFFERENT STARTING POINT (1st PARAMETER) AND OPTIONALLY DIFFERENT DIRECTION (2nd PARAMETER)

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"ground");
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body->frame.t=VECTOR_3D<T>(0,0,0);
    rigid_body->velocity=VECTOR_3D<T>(0,0,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;

    RIGID_BODY_LIST_3D<T>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    for(int i=1;i<=rigid_body_list.Number_Of_Elements();i++) if(!rigid_body_particles.Rigid_Body(i).is_static)
        rigid_body_particles.Rigid_Body(i).Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);


    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Make_Lathe_Chain
//#####################################################################
template<class T,class RW> void SNAKE_EXAMPLE<T,RW>::
Make_Snake_Chain(VECTOR_3D<T> start_point,VECTOR_3D<T> direction,VECTOR_3D<T> shift,QUATERNION<T> orient,int& num_joints,int& num_bodies,std::string name)
{
    RIGID_BODY<TV> *rigid_body = 0;
    int id=0;
    T scale_factor=1;
    T cor=(T).3;
    T cof=(T).6;

    direction.Normalize();
    T link_length=8; // cylinder length
    VECTOR_3D<T> current_position=start_point;
    QUATERNION<T> current_orientation=QUATERNION<T>::Rotation_Quaternion(VECTOR_3D<T>(0,0,1),direction);// rotate the up vector to the direction vector
    // CHANGE NUMBER OF LINKS HERE TO MAKE LONGER
    int number_of_links=8;

    //1
    for(int i=1;i<=number_of_links;i++){
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"Rings_Test/cylinder_revolve",scale_factor);
        rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
        rigid_body->frame.t=current_position+direction*((T).5*link_length);
        rigid_body->frame.r=orient*QUATERNION<T>::Rotation_Quaternion(VECTOR_3D<T>(0,0,1),direction);
        rigid_body->Set_Coefficient_Of_Restitution(cor);
        rigid_body->Set_Coefficient_Of_Friction(cof);
        // HERE'S THE MASS . . . UM, YEAH, GUESS THAT WAS KINDA OBVIOUS
        rigid_body->Set_Mass(rigid_body->mass*10*i); // NOTE: was rigid_body->mass*=10*i which didn't update inertia tensor
        current_position+=direction*link_length;
        rigid_body->name=name;
    }

    ANGLE_JOINT<TV>* joint;
    if(!this->restart){
        for(int i=1;i<=number_of_links;i++) {
            if(i==1) {bend_function_joint=new BEND_FUNCTION_JOINT<TV>();joint=bend_function_joint;}
            else joint=new ANGLE_JOINT<TV>();
            arb->joint_mesh.Add_Joint(joint);
            joint->Set_Joint_To_Child_Frame(FRAME_3D<T>(VECTOR_3D<T>(0,(T).5*link_length,0),QUATERNION<T>((T)(.5*pi),VECTOR_3D<T>(0,1,0))));
            joint->Set_Joint_To_Parent_Frame(FRAME_3D<T>(VECTOR_3D<T>(0,-(T).5*link_length,0),QUATERNION<T>((T)(.5*pi),VECTOR_3D<T>(0,1,0))));
            arb->Add_Articulation(num_bodies+i,num_bodies+i+1,num_joints+i);
        }
        bend_function_joint->Set_Joint_To_Parent_Frame(arb->rigid_bodies_list.rigid_bodies(1)->frame.Inverse()*arb->rigid_bodies_list.rigid_bodies(2)->frame*bend_function_joint->F_cj());


        num_joints+=number_of_links;}
    num_bodies+=number_of_links;
}
//#####################################################################
// Function Run
//#####################################################################
template<class T,class RW> void SNAKE_EXAMPLE<T,RW>::
Postprocess_Frame(const int frame)
{
}
//#####################################################################
// Function )Preprocess
//#####################################################################
template<class T,class RW> void SNAKE_EXAMPLE<T,RW>::
Preprocess_Frame(const int frame)
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Postprocess_Frame(frame);
    std::cout<<"RUNNING preprocess . . . .\n";

    bend_function_joint->Evaluate_Function(frame);
/*
    VECTOR_3D<T> angular_momentum_new=VECTOR_3D<T>(0,current_frame*.1,0);
    arb->rigid_bodies_list.rigid_bodies(2)->Angular_Momentum()=angular_momentum_new;

    if(flip_frame) current_frame++;else current_frame--;
    if(current_frame==10 || current_frame==-10) flip_frame=!flip_frame;
*/
}

namespace PhysBAM
{
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    template class SNAKE_EXAMPLE<double,float>;
#endif
    template class SNAKE_EXAMPLE<float,float>;
}
