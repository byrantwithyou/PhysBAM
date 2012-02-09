//#####################################################################
// Copyright 2004, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "JOINT_TYPES_EXAMPLE.h"

using namespace PhysBAM;

template<class T>
JOINT_TYPES_EXAMPLE<T>::JOINT_TYPES_EXAMPLE(int parameter) : ARTICULATED_RIGID_BODIES_3D_EXAMPLE<T,T>()
{
    last_frame=500;
    frame_rate=24;
    output_directory="Joint_Types/output";
    std::cout << "Frame rate: "<<frame_rate<<std::endl;
}

template<class T>
JOINT_TYPES_EXAMPLE<T>::~JOINT_TYPES_EXAMPLE()
{}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void JOINT_TYPES_EXAMPLE<T>::
Initialize_Rigid_Bodies()
{
    Fixed_Joint(-12);
    Ball_Joint(-6);
    Pin_Joint(0);
    Hinge_Joint(6);
    Universal_Joint(12);
    Gimbal_Joint(18);

    RIGID_BODY<TV> *rigid_body = 0;
    rigid_body=Initialize_Rigid_Body("ground",.3);
    rigid_body->position=VECTOR_3D<T>(0,-3,0);
    rigid_body->velocity=VECTOR_3D<T>(0,0,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
}
//#####################################################################
// Initialization Functions
//#####################################################################
template<class T> void JOINT_TYPES_EXAMPLE<T>::
Fixed_Joint(const T shift)
{
    arb.joint_mesh.joint_description_list.Add_Joint_Description();
    JOINT<TV>* joint=arb.joint_mesh.joints(arb.joint_mesh.Add_Joint());

    joint->joint_description.prismatic_component=new JOINT_PRISMATIC_COMPONENT_3D<T>(0);
    joint->joint_description.prismatic_component->translation_0d=VECTOR_3D<T>(0,0,0);

    // constrain everything
    joint->joint_description.constraints=new JOINT_CONSTRAINTS_3D<T>();
    joint->joint_description.constraints->Set_Twist_Constraints(true,0,0);
    joint->joint_description.constraints->Set_Rotation_Constraints_Max_Min(true,0,0,0,0);

    RIGID_BODY<TV> *rigid_body = 0;
    rigid_body=Initialize_Rigid_Body("subdivided_box");
    rigid_body->position=VECTOR_3D<T>(shift,2,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->Set_Name("parent");
    joint->joint_description.Set_Joint_To_Parent_Frame(FRAME<T>(VECTOR_3D<T>(1,1,1)));
    
    rigid_body=Initialize_Rigid_Body("subdivided_box");
    rigid_body->position=VECTOR_3D<T>(shift+2,4,2);
    rigid_body->velocity=VECTOR_3D<T>(0,0,0);
    rigid_body->Angular_Momentum()=VECTOR_3D<T>(5,-3,10);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->Set_Name("child");
    joint->joint_description.Set_Joint_To_Child_Frame(FRAME<T>(VECTOR_3D<T>(-1,-1,-1)));

    arb.Add_Articulation(1,2,1);
}
template<class T> void JOINT_TYPES_EXAMPLE<T>::
Ball_Joint(const T shift)
{
    arb.joint_mesh.joint_description_list.Add_Joint_Description();
    JOINT<TV>* joint=arb.joint_mesh.joints(arb.joint_mesh.Add_Joint());

    joint->joint_description.prismatic_component=new JOINT_PRISMATIC_COMPONENT_3D<T>(0);
    joint->joint_description.prismatic_component->translation_0d=VECTOR_3D<T>(0,0,0);

    // no constraints

    RIGID_BODY<TV> *rigid_body = 0;

    rigid_body=Initialize_Rigid_Body("subdivided_box");
    rigid_body->position=VECTOR_3D<T>(shift,2,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->Set_Name("parent");
    joint->joint_description.Set_Joint_To_Parent_Frame(FRAME<T>(VECTOR_3D<T>(1,1,1)));
    
    rigid_body=Initialize_Rigid_Body("subdivided_box");
    rigid_body->position=VECTOR_3D<T>(shift+2,4,2);
    rigid_body->velocity=VECTOR_3D<T>(0,0,0);
    rigid_body->Angular_Momentum()=VECTOR_3D<T>(5,0,10);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->Set_Name("child");
    joint->joint_description.Set_Joint_To_Child_Frame(FRAME<T>(VECTOR_3D<T>(-1,-1,-1)));

    arb.Add_Articulation(1,2,1);
}
template<class T> void JOINT_TYPES_EXAMPLE<T>::
Pin_Joint(const T shift)
{
    arb.joint_mesh.joint_description_list.Add_Joint_Description();
    JOINT<TV>* joint=arb.joint_mesh.joints(arb.joint_mesh.Add_Joint());

    joint->joint_description.prismatic_component=new JOINT_PRISMATIC_COMPONENT_3D<T>(0);
    joint->joint_description.prismatic_component->translation_0d=VECTOR_3D<T>(0.1,0,0);

    // constrain theta and phi, not twist
    joint->joint_description.constraints=new JOINT_CONSTRAINTS_3D<T>();
    joint->joint_description.constraints->Set_Rotation_Constraints_Max_Min(true,0,0,0,0);

    RIGID_BODY<TV> *rigid_body = 0;

    rigid_body=Initialize_Rigid_Body("subdivided_box");
    rigid_body->position=VECTOR_3D<T>(shift,2,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->Set_Name("parent");
    joint->joint_description.Set_Joint_To_Parent_Frame(FRAME<T>(VECTOR_3D<T>(1,0,0)));
    
    rigid_body=Initialize_Rigid_Body("subdivided_box");
    rigid_body->position=VECTOR_3D<T>(shift+2.1,2,0);
    rigid_body->velocity=VECTOR_3D<T>(0,0,0);
    rigid_body->Angular_Momentum()=VECTOR_3D<T>(10,0,10);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->Set_Name("child");
    joint->joint_description.Set_Joint_To_Child_Frame(FRAME<T>(VECTOR_3D<T>(-1,0,0)));

    arb.Add_Articulation(1,2,1);
}

template<class T> void JOINT_TYPES_EXAMPLE<T>::
Hinge_Joint(const T shift)
{
    arb.joint_mesh.joint_description_list.Add_Joint_Description();
    JOINT<TV>* joint=arb.joint_mesh.joints(arb.joint_mesh.Add_Joint());

    joint->joint_description.prismatic_component=new JOINT_PRISMATIC_COMPONENT_3D<T>(0);
    joint->joint_description.prismatic_component->translation_0d=VECTOR_3D<T>(0,0,0);
    joint->joint_description.constraints=new JOINT_CONSTRAINTS_3D<T>();
    joint->joint_description.constraints->Set_Twist_Constraints(true,0,0);
    joint->joint_description.constraints->Set_Rotation_Constraints_Max_Min(true,-pi,pi,0,0);

    arb.joint_mesh.joint_description_list.Add_Joint_Description();
    JOINT<TV>* joint2=arb.joint_mesh.joints(arb.joint_mesh.Add_Joint());
    
    RIGID_BODY<TV> *rigid_body = 0;

    rigid_body=Initialize_Rigid_Body("subdivided_box");
    rigid_body->position=VECTOR_3D<T>(shift,2,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->Set_Name("parent");
    joint->joint_description.Set_Joint_To_Parent_Frame(FRAME<T>(VECTOR_3D<T>(1,1,1)));
    joint2->joint_description.Set_Joint_To_Parent_Frame(FRAME<T>(VECTOR_3D<T>(1,1,-1)));
    
    rigid_body=Initialize_Rigid_Body("subdivided_box");
    rigid_body->position=VECTOR_3D<T>(shift+2,4,0);
    rigid_body->velocity=VECTOR_3D<T>(0,0,0);
    rigid_body->Angular_Momentum()=VECTOR_3D<T>(0,0,10);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->Set_Name("child");
    joint->joint_description.Set_Joint_To_Child_Frame(FRAME<T>(VECTOR_3D<T>(-1,-1,1)));
    joint2->joint_description.Set_Joint_To_Child_Frame(FRAME<T>(VECTOR_3D<T>(-1,-1,-1)));

    arb.Add_Articulation(1,2,1);
    arb.Add_Articulation(1,2,2);
}

template<class T> void JOINT_TYPES_EXAMPLE<T>::
Universal_Joint(const T shift)
{
    arb.joint_mesh.joint_description_list.Add_Joint_Description();
    JOINT<TV>* joint=arb.joint_mesh.joints(arb.joint_mesh.Add_Joint());

    joint->joint_description.prismatic_component=new JOINT_PRISMATIC_COMPONENT_3D<T>(0);
    joint->joint_description.prismatic_component->translation_0d=VECTOR_3D<T>(0,0,0);

    // constrain twist
    joint->joint_description.constraints=new JOINT_CONSTRAINTS_3D<T>();
    joint->joint_description.constraints->Set_Twist_Constraints(true,0,0);
    joint->joint_description.constraints->Set_Rotation_Constraints_Max_Min(true,-pi,pi,-pi,pi);

    RIGID_BODY<TV> *rigid_body = 0;
    rigid_body=Initialize_Rigid_Body("subdivided_box");
    rigid_body->position=VECTOR_3D<T>(shift,2,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->Set_Name("parent");
    joint->joint_description.Set_Joint_To_Parent_Frame(FRAME<T>(VECTOR_3D<T>(1,1,1),QUATERNION<T>(-pi/2,VECTOR_3D<T>(0,1,0))));
    
    rigid_body=Initialize_Rigid_Body("subdivided_box");
    rigid_body->position=VECTOR_3D<T>(shift+2,4,2);
    rigid_body->velocity=VECTOR_3D<T>(0,0,0);
    rigid_body->Angular_Momentum()=VECTOR_3D<T>(5,-3,10);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->Set_Name("child");
    joint->joint_description.Set_Joint_To_Child_Frame(FRAME<T>(VECTOR_3D<T>(-1,-1,-1),QUATERNION<T>(-pi/2,VECTOR_3D<T>(0,1,0))));

    arb.Add_Articulation(1,2,1);
}

template<class T> void JOINT_TYPES_EXAMPLE<T>::
Gimbal_Joint(const T shift)
{
    arb.joint_mesh.joint_description_list.Add_Joint_Description();
    JOINT<TV>* joint=arb.joint_mesh.joints(arb.joint_mesh.Add_Joint());
    joint->joint_description.prismatic_component=new JOINT_PRISMATIC_COMPONENT_3D<T>(0);
    joint->joint_description.prismatic_component->translation_0d=VECTOR_3D<T>(0,0,0);

    arb.joint_mesh.joint_description_list.Add_Joint_Description();
    JOINT<TV>* joint2=arb.joint_mesh.joints(arb.joint_mesh.Add_Joint());

    arb.joint_mesh.joint_description_list.Add_Joint_Description();
    JOINT<TV>* joint3=arb.joint_mesh.joints(arb.joint_mesh.Add_Joint());

    // constrain all but twist
    joint->joint_description.constraints=new JOINT_CONSTRAINTS_3D<T>();
    joint->joint_description.constraints->Set_Twist_Constraints(false);
    joint->joint_description.constraints->Set_Rotation_Constraints_Max_Min(true,0,0,0,0);
    // constrain all but twist
    joint2->joint_description.constraints=new JOINT_CONSTRAINTS_3D<T>();
    joint2->joint_description.constraints->Set_Twist_Constraints(false);
    joint2->joint_description.constraints->Set_Rotation_Constraints_Max_Min(true,0,0,0,0);
    // constrain all but twist
    joint3->joint_description.constraints=new JOINT_CONSTRAINTS_3D<T>();
    joint3->joint_description.constraints->Set_Twist_Constraints(false);
    joint3->joint_description.constraints->Set_Rotation_Constraints_Max_Min(true,0,0,0,0);

    RIGID_BODY<TV> *rigid_body = 0;
    rigid_body=Initialize_Rigid_Body("subdivided_box");
    rigid_body->position=VECTOR_3D<T>(shift,2,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->Set_Name("parent");
    joint->joint_description.Set_Joint_To_Parent_Frame(FRAME<T>(VECTOR_3D<T>(1,1,1)));
    joint2->joint_description.Set_Joint_To_Parent_Frame(FRAME<T>(VECTOR_3D<T>(1,1,1),QUATERNION<T>(-pi/2,VECTOR_3D<T>(0,1,0))));
    joint3->joint_description.Set_Joint_To_Parent_Frame(FRAME<T>(VECTOR_3D<T>(1,1,1),QUATERNION<T>(-pi/2,VECTOR_3D<T>(1,0,0))));
    
    rigid_body=Initialize_Rigid_Body("subdivided_box");
    rigid_body->position=VECTOR_3D<T>(shift+2,4,2);
    rigid_body->velocity=VECTOR_3D<T>(0,0,0);
    rigid_body->Angular_Momentum()=VECTOR_3D<T>(5,0,10);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    rigid_body->Set_Name("child");
    joint->joint_description.Set_Joint_To_Child_Frame(FRAME<T>(VECTOR_3D<T>(-1,-1,-1)));
    joint2->joint_description.Set_Joint_To_Child_Frame(FRAME<T>(VECTOR_3D<T>(-1,-1,-1),QUATERNION<T>(-pi/2,VECTOR_3D<T>(0,1,0))));
    joint3->joint_description.Set_Joint_To_Child_Frame(FRAME<T>(VECTOR_3D<T>(-1,-1,-1),QUATERNION<T>(-pi/2,VECTOR_3D<T>(1,0,0))));

    arb.Add_Articulation(0,1,0);
    arb.Add_Articulation(0,1,1);
    arb.Add_Articulation(0,1,2);
}

namespace PhysBAM
{
    template JOINT_TYPES_EXAMPLE<double>;
    template JOINT_TYPES_EXAMPLE<float>;
}
