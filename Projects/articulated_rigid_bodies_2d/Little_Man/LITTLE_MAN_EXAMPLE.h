//#####################################################################
// Copyright 2004, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LITTLE_MAN_EXAMPLE
//##################################################################### 
#ifndef __LITTLE_MAN_EXAMPLE__
#define __LITTLE_MAN_EXAMPLE__

#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T,class RW>
class LITTLE_MAN_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;
    using BASE::restart;using BASE::restart_frame;using BASE::output_directory;
    using BASE::solids_parameters;using BASE::write_last_frame;using BASE::data_directory;

    ARTICULATED_RIGID_BODY<TV>* arb;

    LITTLE_MAN_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>(FLUIDS_PARAMETERS_2D<T>::NONE)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=.1;
        //solids_parameters.gravity=0;

        last_frame=3000;
        frame_rate=24;
        output_directory="Little_Man/output";
        std::cout << "Frame rate: "<<frame_rate<<std::endl;

        arb=new ARTICULATED_RIGID_BODY<TV>(this->solids_parameters.rigid_body_parameters.list);
        this->solids_parameters.rigid_body_parameters.Set_Articulated_Rigid_Body(arb);
        arb->Set_Do_Final_Pass(false);
    }

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    RIGID_BODY<TV>* rigid_body=0;
    Make_A_Body();

    std::cout <<"start initing bodies\n";
    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.1);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body->frame.t=VECTOR_2D<T>(0,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(1);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;

    RIGID_BODY_LIST_2D<T>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    for(int i=0;i<rigid_body_list.Number_Of_Elements();i++) if(!rigid_body_particles.Rigid_Body(i).is_static)
        rigid_body_particles.Rigid_Body(i).Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,(T)0);

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::Initialize_Bodies();
    std::cout <<"done initing bodies\n";
}
//#####################################################################
// Function Make_A_Body
//#####################################################################
void Make_A_Body()
{
    T friction=1,epsilon=0;
    T k_p=10;
//    FRAME_2D<T> transform(VECTOR_2D<T>(0,2),COMPLEX<T>::Unit_Polar(0.1));
    FRAME_2D<T> transform(VECTOR_2D<T>(0,2),COMPLEX<T>::Unit_Polar(0));
    RIGID_BODY<TV> *rigid_body = 0;

//1
    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body->frame.t=VECTOR_2D<T>(1,1);
    rigid_body->frame.r=COMPLEX<T>::Unit_Polar(pi/2);
    //rigid_body->frame.r=COMPLEX<T>::Unit_Polar(-pi/4);
    //rigid_body->Angular_Momentum()=-10;
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Name("calfright");
    
//2 
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    //rigid_body->frame.t=VECTOR_2D<T>(sqrt(2.)+1,1);
    rigid_body->frame.r=COMPLEX<T>::Unit_Polar(pi/2);
    rigid_body->frame.t=VECTOR_2D<T>(1,3);
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Name("thighright");

//3
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body->frame.t=VECTOR_2D<T>(-1,1);
    rigid_body->frame.r=COMPLEX<T>::Unit_Polar(pi/2);
    //rigid_body->frame.r=COMPLEX<T>::Unit_Polar(-pi/4);
    //rigid_body->Angular_Momentum()=-10;
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Name("calfleft");
   
//4 
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    //rigid_body->frame.t=VECTOR_2D<T>(sqrt(2.)+1,1);
    rigid_body->frame.r=COMPLEX<T>::Unit_Polar(pi/2);
    rigid_body->frame.t=VECTOR_2D<T>(-1,3);
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Name("thighleft");
   
//5 
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    //rigid_body->frame.t=VECTOR_2D<T>(sqrt(2.)+1,1);
    rigid_body->frame.t=VECTOR_2D<T>(0,4);
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Name("hips");

//6
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01*.75);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body->frame.t=VECTOR_2D<T>(0,4.75);
    rigid_body->frame.r=COMPLEX<T>::Unit_Polar(pi/2);
    //rigid_body->frame.r=COMPLEX<T>::Unit_Polar(-pi/4);
    //rigid_body->Angular_Momentum()=-10;
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Name("torsobottom");
   
//7 
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01*.75);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    //rigid_body->frame.t=VECTOR_2D<T>(sqrt(2.)+1,1);
    rigid_body->frame.r=COMPLEX<T>::Unit_Polar(pi/2);
    rigid_body->frame.t=VECTOR_2D<T>(0,6.25);
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Name("torsotop");
   
//8 
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01*1.25);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    //rigid_body->frame.t=VECTOR_2D<T>(sqrt(2.)+1,1);
    rigid_body->frame.t=VECTOR_2D<T>(0,7);
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Name("shoulders");

//9
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body->frame.t=VECTOR_2D<T>(1.25,4);
    rigid_body->frame.r=COMPLEX<T>::Unit_Polar(pi/2);
    //rigid_body->frame.r=COMPLEX<T>::Unit_Polar(-pi/4);
    //rigid_body->Angular_Momentum()=-10;
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Name("forearmright");
    
//10 
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    //rigid_body->frame.t=VECTOR_2D<T>(sqrt(2.)+1,1);
    rigid_body->frame.r=COMPLEX<T>::Unit_Polar(pi/2);
    rigid_body->frame.t=VECTOR_2D<T>(1.25,6);
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Name("bicepright");

//11
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body->frame.t=VECTOR_2D<T>(-1.25,4);
    rigid_body->frame.r=COMPLEX<T>::Unit_Polar(pi/2);
    //rigid_body->frame.r=COMPLEX<T>::Unit_Polar(-pi/4);
    //rigid_body->Angular_Momentum()=-10;
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Name("forearmleft");
   
//12
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    //rigid_body->frame.t=VECTOR_2D<T>(sqrt(2.)+1,1);
    rigid_body->frame.r=COMPLEX<T>::Unit_Polar(pi/2);
    rigid_body->frame.t=VECTOR_2D<T>(-1.25,6);
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Name("bicepleft");
    
//13
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.005);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    //rigid_body->frame.t=VECTOR_2D<T>(sqrt(2.)+1,1);
    rigid_body->frame.r=COMPLEX<T>::Unit_Polar(pi/2);
    rigid_body->frame.t=VECTOR_2D<T>(0,7.5);
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Name("neck");
   
//14
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/circle",.5);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    //rigid_body->frame.t=VECTOR_2D<T>(sqrt(2.)+1,1);
    rigid_body->frame.t=VECTOR_2D<T>(0,8.5);
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(friction);
    rigid_body->Set_Mass(rigid_body->mass*.25);
    rigid_body->Set_Name("head");

    // overall translation and rotation
    for(int i=0;i<solids_parameters.rigid_body_parameters.list.Number_Of_Elements();i++)
        solids_parameters.rigid_body_parameters.list(i)->frame=transform*solids_parameters.rigid_body_parameters.list(i)->frame;

    JOINT<TV>* joint;JOINT_FUNCTION<TV>* jfunc;

//right knee
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(-1,0)));
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(1,0)));
    jfunc=new JOINT_FUNCTION<TV>(joint,arb->rigid_bodies_list.rigid_bodies(2),arb->rigid_bodies_list.rigid_bodies(1));jfunc->Set_Desired_Angle(jfunc->Angle());
    jfunc->Set_k_p(k_p);
    joint->Set_Joint_Function(jfunc);
    id=arb->joint_mesh.Add_Joint(joint);arb->Add_Articulation(2,1,id);

//left knee
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(-1,0)));
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(1,0)));
    jfunc=new JOINT_FUNCTION<TV>(joint,arb->rigid_bodies_list.rigid_bodies(4),arb->rigid_bodies_list.rigid_bodies(3));jfunc->Set_Desired_Angle(jfunc->Angle());
    jfunc->Set_k_p(k_p);
    joint->Set_Joint_Function(jfunc);
    id=arb->joint_mesh.Add_Joint(joint);arb->Add_Articulation(4,3,id);

//right hip
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(1,0)));
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(1,0)));
    jfunc=new JOINT_FUNCTION<TV>(joint,arb->rigid_bodies_list.rigid_bodies(5),arb->rigid_bodies_list.rigid_bodies(2));jfunc->Set_Desired_Angle(jfunc->Angle()-0.3);
    jfunc->Set_k_p(k_p);
    joint->Set_Joint_Function(jfunc);
    id=arb->joint_mesh.Add_Joint(joint);arb->Add_Articulation(5,2,id);

//left hip
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(-1,0)));
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(1,0)));
    jfunc=new JOINT_FUNCTION<TV>(joint,arb->rigid_bodies_list.rigid_bodies(5),arb->rigid_bodies_list.rigid_bodies(4));jfunc->Set_Desired_Angle(jfunc->Angle()+0.3);
    jfunc->Set_k_p(k_p);
    joint->Set_Joint_Function(jfunc);
    id=arb->joint_mesh.Add_Joint(joint);arb->Add_Articulation(5,4,id);

//pelvis
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(-.75,0)));
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(0,0)));
    jfunc=new JOINT_FUNCTION<TV>(joint,arb->rigid_bodies_list.rigid_bodies(6),arb->rigid_bodies_list.rigid_bodies(5));jfunc->Set_Desired_Angle(jfunc->Angle());
    jfunc->Set_k_p(k_p);
    joint->Set_Joint_Function(jfunc);
    id=arb->joint_mesh.Add_Joint(joint);arb->Add_Articulation(6,5,id);

//stomach
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(.75,0)));
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(-.75,0)));
    jfunc=new JOINT_FUNCTION<TV>(joint,arb->rigid_bodies_list.rigid_bodies(6),arb->rigid_bodies_list.rigid_bodies(7));jfunc->Set_Desired_Angle(jfunc->Angle());
    jfunc->Set_k_p(k_p);
    joint->Set_Joint_Function(jfunc);
    id=arb->joint_mesh.Add_Joint(joint);arb->Add_Articulation(6,7,id);

//sternum
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(.75,0)));
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(0,0)));
    jfunc=new JOINT_FUNCTION<TV>(joint,arb->rigid_bodies_list.rigid_bodies(7),arb->rigid_bodies_list.rigid_bodies(8));jfunc->Set_Desired_Angle(jfunc->Angle());
    jfunc->Set_k_p(k_p);
    joint->Set_Joint_Function(jfunc);
    id=arb->joint_mesh.Add_Joint(joint);arb->Add_Articulation(7,8,id);

//right elbow
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(-1,0)));
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(1,0)));
    jfunc=new JOINT_FUNCTION<TV>(joint,arb->rigid_bodies_list.rigid_bodies(10),arb->rigid_bodies_list.rigid_bodies(9));jfunc->Set_Desired_Angle(jfunc->Angle());
    jfunc->Set_k_p(k_p);
    joint->Set_Joint_Function(jfunc);
    id=arb->joint_mesh.Add_Joint(joint);arb->Add_Articulation(10,9,id);

//left elbow
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(-1,0)));
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(1,0)));
    jfunc=new JOINT_FUNCTION<TV>(joint,arb->rigid_bodies_list.rigid_bodies(12),arb->rigid_bodies_list.rigid_bodies(11));jfunc->Set_Desired_Angle(jfunc->Angle());
    jfunc->Set_k_p(k_p);
    joint->Set_Joint_Function(jfunc);
    id=arb->joint_mesh.Add_Joint(joint);arb->Add_Articulation(12,11,id);

//right shoulder
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(1.25,0)));
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(1,0)));
    jfunc=new JOINT_FUNCTION<TV>(joint,arb->rigid_bodies_list.rigid_bodies(8),arb->rigid_bodies_list.rigid_bodies(10));jfunc->Set_Desired_Angle(jfunc->Angle());
    jfunc->Set_k_p(k_p);
    joint->Set_Joint_Function(jfunc);
    id=arb->joint_mesh.Add_Joint(joint);arb->Add_Articulation(8,10,id);

//left shoulder
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(-1.25,0)));
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(1,0)));
    jfunc=new JOINT_FUNCTION<TV>(joint,arb->rigid_bodies_list.rigid_bodies(8),arb->rigid_bodies_list.rigid_bodies(12));jfunc->Set_Desired_Angle(jfunc->Angle());
    jfunc->Set_k_p(k_p);
    joint->Set_Joint_Function(jfunc);
    id=arb->joint_mesh.Add_Joint(joint);arb->Add_Articulation(8,12,id);

#if 0
//neck
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(0,0)));
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(-.5,0)));
    jfunc=new JOINT_FUNCTION<TV>(joint,arb->rigid_bodies_list.rigid_bodies(8),arb->rigid_bodies_list.rigid_bodies(13));jfunc->Set_Desired_Angle(jfunc->Angle());
    jfunc->Set_k_p(k_p);
    joint->Set_Joint_Function(jfunc);
    id=arb->joint_mesh.Add_Joint(joint);arb->Add_Articulation(8,13,id);

//necktop
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(.5,0)));
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(0,-.5)));
    jfunc=new JOINT_FUNCTION<TV>(joint,arb->rigid_bodies_list.rigid_bodies(13),arb->rigid_bodies_list.rigid_bodies(14));jfunc->Set_Desired_Angle(jfunc->Angle());
    jfunc->Set_k_p(k_p);
    joint->Set_Joint_Function(jfunc);
    id=arb->joint_mesh.Add_Joint(joint);arb->Add_Articulation(13,14,id);
#endif
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE{
    // objects only collide with ground
    RIGID_BODY_LIST_2D<T>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    solids_evolution->rigid_body_collisions->collision_manager.Use_Collision_Matrix();
    for(int i=1;i<=rigid_body_list.Number_Of_Elements()-1;i++) for(int j=1;j<=rigid_body_list.Number_Of_Elements()-1;j++)
        solids_evolution->rigid_body_collisions->collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(i,j,false);
}
//#####################################################################
};
}
#endif
