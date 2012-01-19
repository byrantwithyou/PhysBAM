//#####################################################################
// Copyright 2004, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HORSE_EXAMPLE
//##################################################################### 
#ifndef __HORSE_EXAMPLE__
#define __HORSE_EXAMPLE__

#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/RIGID_JOINT.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T,class RW>
class HORSE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<RW,GRID<TV> >
{
    typedef VECTOR<T,3> TV;
public:
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<RW,GRID<TV> >::first_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<RW,GRID<TV> >::last_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<RW,GRID<TV> >::frame_rate;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<RW,GRID<TV> >::restart;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<RW,GRID<TV> >::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<RW,GRID<TV> >::output_directory;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<RW,GRID<TV> >::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<RW,GRID<TV> >::write_last_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<RW,GRID<TV> >::data_directory;
    
    ARTICULATED_RIGID_BODY<TV>* arb;
    RIGID_BODY<TV> *shelf11,*shelf12,*shelf21,*shelf22;
    int current_frame,start_move,end_move;
    T increment;
    int selection;
    int id1,id2;
    bool flip;
    JOINT<TV>* joint[4];

    HORSE_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<RW,GRID<TV> >(FLUIDS_PARAMETERS_3D<T>::NONE)
    {
        //solids_parameters.perform_self_collision=false;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        //solids_parameters.gravity=0;
        //solids_parameters.perform_collision_body_collisions=false;
    //    restart=true;
        restart_frame=30;
        last_frame=2000;
        frame_rate=24*2;
        output_directory="Horse/output";
        std::cout << "Frame rate: "<<frame_rate<<std::endl;

        arb=new ARTICULATED_RIGID_BODY<TV>(this->solids_parameters.rigid_body_parameters.list);
        this->solids_parameters.rigid_body_parameters.Set_Articulated_Rigid_Body(arb);
        
        arb->Set_Iterative_Tolerance((T)1e-4);
        current_frame=0;
        increment=(T).05;
        start_move=5;end_move=40;
        shelf11=shelf12=shelf21=shelf22=0;
    //    arb->Set_Extra_Iterations_Per_Contact_Level_Factor(100);
    //    arb->Set_Extra_Iterations_Per_Shock_Propagation_Level_Factor(100);
    //    arb->Set_Poststabilization_Iterations(100);
        arb->Set_Use_Shock_Propagation(false);
        arb->Set_Do_Final_Pass(false);
        write_last_frame=true;
        selection=1;
        flip=true;
    }

    ~HORSE_EXAMPLE()
    {
        delete arb;
    }

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    int id=0;
    RIGID_BODY<TV> *child_body=0,*rigid_body=0,*parent_body=0;
    //JOINT<TV>* joint[4];
    if(!this->restart){ 
        for(int i=0;i<4;i++){joint[i]=new ANGLE_JOINT<TV>();arb->joint_mesh.Add_Joint(joint[i]);}
    }

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"plank");
    parent_body=arb->rigid_bodies_list.rigid_bodies(id);
    parent_body->frame.t=TV(0,2,0);
    parent_body->Set_Coefficient_Of_Restitution(0.5);
    parent_body->Set_Coefficient_Of_Friction(.5);
    parent_body->Set_Name("body");
    
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"Rings_Test/cylinder_revolve",.5);
    child_body=arb->rigid_bodies_list.rigid_bodies(id);
    child_body->frame.t=TV(-1.25,0,4);
    child_body->Set_Coefficient_Of_Restitution(0.5);
    child_body->Set_Coefficient_Of_Friction(.5);
    child_body->Set_Name("leg");
    joint[0]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(.25,2,0)));
    joint[0]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,0,4)));
    arb->joint_mesh.Add_Articulation(1,2,1);
    joint[0]->Set_Joint_Function(new JOINT_FUNCTION<TV>(joint[0],parent_body,child_body));joint[0]->joint_function->Set_k_p(30);
    joint[0]->joint_function->Set_Desired_Angle(QUATERNION<T>(-.3,TV(1,0,0)));

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"Rings_Test/cylinder_revolve",.5);
    child_body=arb->rigid_bodies_list.rigid_bodies(id);
    child_body->frame.t=TV(-1.25,0,-4);
    child_body->Set_Coefficient_Of_Restitution(0.5);
    child_body->Set_Coefficient_Of_Friction(.5);
    child_body->Set_Name("leg");
    joint[1]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(.25,2,0)));
    joint[1]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,0,-4)));
    if(!this->restart) arb->joint_mesh.Add_Articulation(1,3,2);
    joint[1]->Set_Joint_Function(new JOINT_FUNCTION<TV>(joint[1],parent_body,child_body));joint[1]->joint_function->Set_k_p(30);
    joint[1]->joint_function->Set_Desired_Angle(QUATERNION<T>(.3,TV(1,0,0)));
    
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"Rings_Test/cylinder_revolve",.5);
    child_body=arb->rigid_bodies_list.rigid_bodies(id);
    child_body->frame.t=TV(1.25,0,-4);
    child_body->Set_Coefficient_Of_Restitution(0.5);
    child_body->Set_Coefficient_Of_Friction(.5);
    child_body->Set_Name("leg");
    joint[2]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-.25,2,0)));
    joint[2]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,0,-4)));
    if(!this->restart) arb->joint_mesh.Add_Articulation(1,4,3);
    joint[2]->Set_Joint_Function(new JOINT_FUNCTION<TV>(joint[2],parent_body,child_body));joint[2]->joint_function->Set_k_p(30);
    joint[2]->joint_function->Set_Desired_Angle(QUATERNION<T>(-.3,TV(1,0,0)));
 
   id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"Rings_Test/cylinder_revolve",.5);
    child_body=arb->rigid_bodies_list.rigid_bodies(id);
    child_body->frame.t=TV(1.25,0,4);
    child_body->Set_Coefficient_Of_Restitution(0.5);
    child_body->Set_Coefficient_Of_Friction(.5);
    child_body->Set_Name("leg");
    joint[3]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-.25,2,0)));
    joint[3]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,0,4)));
    if(!this->restart) arb->joint_mesh.Add_Articulation(1,5,4);
    joint[3]->Set_Joint_Function(new JOINT_FUNCTION<TV>(joint[3],parent_body,child_body));joint[3]->joint_function->Set_k_p(30);
    joint[3]->joint_function->Set_Desired_Angle(QUATERNION<T>(.3,TV(1,0,0)));
    
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"ground");
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body->frame.t=TV(0,-2,0);
    rigid_body->velocity=TV(0,0,0);
    rigid_body->Set_Coefficient_Of_Restitution(.5);
    rigid_body->Set_Coefficient_Of_Friction(.5);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;

    RIGID_BODY_LIST_3D<T>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    for(int i=0;i<rigid_body_list.Number_Of_Elements();i++) if(!rigid_body_particles.Rigid_Body(i).is_static)
        rigid_body_particles.Rigid_Body(i).Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<RW,GRID<TV> >::Initialize_Bodies();
    std::cout<<"done initializing example\n";
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(flip){
        joint[0]->joint_function->Set_Desired_Angle(QUATERNION<T>(-.3,TV(1,0,0)));
        joint[1]->joint_function->Set_Desired_Angle(QUATERNION<T>(.3,TV(1,0,0)));
        joint[2]->joint_function->Set_Desired_Angle(QUATERNION<T>(-.3,TV(1,0,0)));
        joint[3]->joint_function->Set_Desired_Angle(QUATERNION<T>(.3,TV(1,0,0)));
    }
    else{
        joint[0]->joint_function->Set_Desired_Angle(QUATERNION<T>(.3,TV(1,0,0)));
        joint[1]->joint_function->Set_Desired_Angle(QUATERNION<T>(-.3,TV(1,0,0)));
        joint[2]->joint_function->Set_Desired_Angle(QUATERNION<T>(.3,TV(1,0,0)));
        joint[3]->joint_function->Set_Desired_Angle(QUATERNION<T>(-.3,TV(1,0,0)));
    }
    if(!(frame%30)) {if(flip)flip=false;else flip=true;}
}
//#####################################################################
};
}
#endif
