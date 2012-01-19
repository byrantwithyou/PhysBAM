//#####################################################################
// Copyright 2004-2005, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLE_MUSCLE_EXAMPLE
//##################################################################### 
#ifndef __SIMPLE_MUSCLE_EXAMPLE__
#define __SIMPLE_MUSCLE_EXAMPLE__

#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_FORCE_CURVE.h>
#include <PhysBAM_Dynamics/Motion/FRAME_TRACK_3D.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <fstream>
#include <Rigid_Bodies/CONSTRAINED_POINT_IN_RIGID_BODY.h>
namespace PhysBAM{

template<class T,class RW>
class SIMPLE_MUSCLE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    typedef CONSTRAINED_POINT_IN_RIGID_BODY<T,VECTOR<T,3> > T_CONSTRAINED_POINT_IN_RIGID_BODY;

    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;
    using BASE::restart;using BASE::restart_frame;using BASE::output_directory;
    using BASE::solids_parameters;using BASE::write_last_frame;using BASE::data_directory;

    ARTICULATED_RIGID_BODY<TV>* arb;
    MUSCLE<T,GRID<TV> >* muscle;
    bool add_ground;
    PARAMETER_LIST parameter_list;
    T peak_force;

    SIMPLE_MUSCLE_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(0,FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >::NONE)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=(T).1;
        //solids_parameters.gravity=0;
        //solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=1;

        last_frame=3000;
        frame_rate=24;
        output_directory="Simple_Muscle/output";
        std::cout << "Frame rate: "<<frame_rate<<std::endl;

        arb=new ARTICULATED_RIGID_BODY<TV>(this->solids_parameters.rigid_body_parameters.list);
        this->solids_parameters.rigid_body_parameters.Set_Articulated_Rigid_Body(arb);

        arb->Set_Do_Final_Pass(false);
        arb->Use_Epsilon_Scale(false);
        arb->Set_Use_Shock_Propagation(false);

        arb->Use_Muscle_Actuators();

        muscle=0;
        
        ARB_PARAMETERS::Read_Common_Parameters("Simple_Muscle/example.param",*this,parameter_list);
        peak_force=parameter_list.Get_Parameter("peak_force",(T)10);
    }

    ~SIMPLE_MUSCLE_EXAMPLE()
    {
        delete arb;
    }

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    arb->muscle_list->muscle_force_curve.Initialize(data_directory); // initialize here rather than constructor since data directory might be set after constructor

    add_ground=true;
    //Arm();
    Simple_Muscle_Across_Joint();
    //Stand_Test();

    if(add_ground){
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground",(T).1);
        RIGID_BODY<TV>* rigid_body=arb->rigid_body_list.rigid_bodies(id);
        rigid_body->frame.t=VECTOR<T,3>(0,0,0);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("ground");
        rigid_body->is_static=true;
        rigid_body->add_to_spatial_partition=false;}

    std::cout <<"\nnumber of muscles is: "<<arb->muscle_list->muscles.m<<std::endl;

    RIGID_BODY_LIST<T,VECTOR<T,3> >& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    for(int i=1;i<=rigid_body_list.Number_Of_Elements();i++) if(!rigid_body_particles.Rigid_Body(i).is_static)
       rigid_body_particles.Rigid_Body(i).Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,(T)0);

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
// Function Simple_Muscle_Across_Joint
//#####################################################################
MUSCLE<T,GRID<TV> >* Add_Basic_Muscle(const std::string& name,RIGID_BODY<TV>& origin_body,const VECTOR<T,3>& origin,RIGID_BODY<TV>& insertion_body,const VECTOR<T,3>& insertion)
{
    MUSCLE<T,GRID<TV> >* muscle=new MUSCLE<T,GRID<TV> >(arb->muscle_list->muscle_force_curve);
    muscle->Set_Name(name);
    muscle->Set_Attachment_Point_1(new T_CONSTRAINED_POINT_IN_RIGID_BODY(origin_body,origin));
    muscle->Set_Attachment_Point_2(new T_CONSTRAINED_POINT_IN_RIGID_BODY(insertion_body,insertion));
    T total_length=muscle->Total_Length();
    muscle->Set_Optimal_Length((T).8*total_length);
    muscle->Set_Tendon_Slack_Length((T).2*total_length);
    muscle->Set_Peak_Force(peak_force);
    muscle->Set_Max_Shortening_Velocity(1);
    arb->muscle_list->Add_Muscle(muscle);
    return muscle;
}

void Simple_Muscle_Across_Joint()
{
    RIGID_BODY_LIST<T,VECTOR<T,3> >& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    assert(rigid_body_list.Number_Of_Elements()==0);

    add_ground=false;
    solids_parameters.gravity=0;

    int num_bodies=parameter_list.Get_Parameter("num_bodies",(int)2);
    bool use_bend_joint=parameter_list.Get_Parameter("use_bend_joint",false);
    VECTOR<T,3> plank_rescale=parameter_list.Get_Parameter("plank_rescale",VECTOR<T,3>(1,1,1));
    T target_angle=parameter_list.Get_Parameter("target_angle",-(T)3*pi/4);
    T initial_angle=parameter_list.Get_Parameter("initial_angle",-(T)pi/2);
    T k_p=parameter_list.Get_Parameter("k_p",(T)100);

    FRAME_TRACK_3D<T>* frame_track=0;
    if(parameter_list.Get_Parameter("use_frame_track",false)){
        int samples=1000;T period=2;
        frame_track=new FRAME_TRACK_3D<T>(samples,0,period);frame_track->periodic=true;
        for(int i=0;i<samples;i++){
            frame_track->trajectory(i)=FRAME_3D<T>(QUATERNION<T>(initial_angle+(target_angle-initial_angle)*0.5*(1-cos(2*pi*(i-1)/(samples-1))),VECTOR<T,3>(1,0,0)));}
        for(int i=0;i<arb->joint_mesh.joints.m;i++) arb->joint_mesh.joints(i)->joint_function->track=frame_track;
    }

    for(int i=0;i<num_bodies;i++){
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/plank",(T).2,true,false,false);
        RIGID_BODY<TV>* rigid_body=arb->rigid_body_list.rigid_bodies(id);
        rigid_body->triangulated_surface->Rescale(plank_rescale.x,plank_rescale.y,plank_rescale.z);
        rigid_body->frame.r=QUATERNION<T>(pi/2,VECTOR<T,3>(0,1,0));
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name(STRING_UTILITIES::string_sprintf("body_%d",i));
        rigid_body->Set_Mass(1);

        if(i>1){
            JOINT<TV>* joint=0;
            if(use_bend_joint) joint=new ANGLE_JOINT<TV>();
            else joint=new POINT_JOINT<TV>();
            joint->Set_Joint_To_Parent_Frame(FRAME_3D<T>(VECTOR<T,3>(0,0,1)));
            joint->Set_Joint_To_Child_Frame(FRAME_3D<T>(VECTOR<T,3>(0,0,-1)));
            int joint_id=arb->joint_mesh.Add_Joint(joint);
            arb->joint_mesh.Add_Articulation(i-1,i,joint_id);
            JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(joint,rigid_body_list(i-1),rigid_body_list(i));
            joint->Set_Joint_Function(jfunc);
            if(frame_track) jfunc->track=frame_track;
            else jfunc->Set_Target_Angle(QUATERNION<T>(target_angle,VECTOR<T,3>(1,0,0)));
            //jfunc->Set_Target_Angle(QUATERNION<T>::From_Euler_Angles((T)-pi/2,0.5,0));
            jfunc->muscle_control=true;
            joint->joint_function->Set_k_p(k_p);
            joint->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(initial_angle,VECTOR<T,3>(1,0,0))));
        }
    }

    rigid_body_particles.Rigid_Body(1).frame.r=QUATERNION<T>(-pi/2,VECTOR<T,3>(0,0,1))*rigid_body_particles.Rigid_Body(1).frame.r;


    for(int i=2;i<=num_bodies;i++){
        std::string suffix=STRING_UTILITIES::string_sprintf("_j%d",i-1);
        RIGID_BODY<TV> *parent=rigid_body_list(i-1),*child=rigid_body_list(i);

        Add_Basic_Muscle("flexor"+suffix,*parent,VECTOR<T,3>(0,0.05,-0.3),*child,VECTOR<T,3>(0,0.05,-0.2));
        if(parameter_list.Get_Parameter("add_extensor",true)){
            MUSCLE<T,GRID<TV> > *extensor=Add_Basic_Muscle("extensor"+suffix,*parent,VECTOR<T,3>(0,-0.05,0.5),*child,VECTOR<T,3>(0,-0.05,-0.5));
            extensor->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*parent,VECTOR<T,3>(0,-0.2,1)));
            extensor->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*child,VECTOR<T,3>(0,-0.2,-1)));}
        if(parameter_list.Get_Parameter("add_redundant_muscles",false)){ // EXTRA MUSCLES
            Add_Basic_Muscle("redundant_flexor"+suffix,*parent,VECTOR<T,3>(0,0.05,0.5),*child,VECTOR<T,3>(0,0.05,0.3));
            MUSCLE<T,GRID<TV> > *redundant_extensor=Add_Basic_Muscle("redundant_extensor"+suffix,*parent,VECTOR<T,3>(0,-0.05,-0.5),*child,VECTOR<T,3>(0,-0.05,0.5));
            redundant_extensor->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*parent,VECTOR<T,3>(0,-0.5,1)));
            redundant_extensor->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*child,VECTOR<T,3>(0,-0.4,-1)));}
        if(parameter_list.Get_Parameter("add_diagonal_muscles",false)){
            MUSCLE<T,GRID<TV> > *diag1=Add_Basic_Muscle("diag1"+suffix,*parent,VECTOR<T,3>(0.2,0.05,-0.5),*child,VECTOR<T,3>(-0.2,0.05,-0.5));
            MUSCLE<T,GRID<TV> > *diag2=Add_Basic_Muscle("diag2"+suffix,*parent,VECTOR<T,3>(-0.2,0.05,-0.5),*child,VECTOR<T,3>(0.2,0.05,-0.5));

            MUSCLE<T,GRID<TV> > *side1=Add_Basic_Muscle("side1"+suffix,*parent,VECTOR<T,3>(0.2,0,-0.5),*child,VECTOR<T,3>(0.2,0,0.7));
    //        side1->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*parent,VECTOR<T,3>(0.4,0,0.95)));
            side1->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*child,VECTOR<T,3>(0.4,0,-0.95)));
            MUSCLE<T,GRID<TV> > *side2=Add_Basic_Muscle("side2"+suffix,*parent,VECTOR<T,3>(-0.2,0,-0.5),*child,VECTOR<T,3>(-0.2,0,0.7));
    //        side2->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*parent,VECTOR<T,3>(-0.4,0,0.95)));
            side2->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*child,VECTOR<T,3>(-0.4,0,-0.95)));
        }
    }

    if(parameter_list.Get_Parameter("add_multijoint_muscle",false)){
        MUSCLE<T,GRID<TV> > *muscle=Add_Basic_Muscle("multijoint",*rigid_body_list(1),VECTOR<T,3>(0,0.05,-0.5),*rigid_body_list(rigid_body_list.rigid_bodies.m),VECTOR<T,3>(0,0.05,-0.5));
    }

    arb->Update_With_Breadth_First_Directed_Graph(1);

    if(parameter_list.Get_Parameter("add_obstacle",false)){
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box",0.3);
        RIGID_BODY<TV>* rigid_body=arb->rigid_body_list.rigid_bodies(id);
        rigid_body->frame.t=VECTOR<T,3>(1,3,0);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("obstacle");
        rigid_body->Set_Mass(1);
        rigid_body->is_static=true;
    }
}
//#####################################################################
// Function Stand_Test
//#####################################################################
void Stand_Test()
{
    RIGID_BODY_LIST<T,VECTOR<T,3> >& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    assert(rigid_body_list.Number_Of_Elements()==0);

    add_ground=parameter_list.Get_Parameter("add_ground",false);
    int num_bodies=parameter_list.Get_Parameter("num_bodies",(int)2);
    T k_p=parameter_list.Get_Parameter("k_p",(T)100);

    for(int i=0;i<num_bodies;i++){
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/plank",(T).2,true,false,false);
        RIGID_BODY<TV>* rigid_body=arb->rigid_body_list.rigid_bodies(id);
        rigid_body->frame.r=QUATERNION<T>(pi/2,VECTOR<T,3>(0,1,0));
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name(STRING_UTILITIES::string_sprintf("body_%d",i));
        rigid_body->Set_Mass(1);

        if(i>1){
            JOINT<TV>* joint=new ANGLE_JOINT<TV>();
            joint->Set_Joint_To_Parent_Frame(FRAME_3D<T>(VECTOR<T,3>(0,0,1)));
            joint->Set_Joint_To_Child_Frame(FRAME_3D<T>(VECTOR<T,3>(0,0,-1)));
            int joint_id=arb->joint_mesh.Add_Joint(joint);
            arb->joint_mesh.Add_Articulation(i-1,i,joint_id);
            JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(joint,rigid_body_list(i-1),rigid_body_list(i));
            joint->Set_Joint_Function(jfunc);
            jfunc->Set_Target_Angle(QUATERNION<T>(-(T)3*pi/4,VECTOR<T,3>(1,0,0)));
            //jfunc->Set_Target_Angle(QUATERNION<T>::From_Euler_Angles((T)-pi/2,0.5,0));
            jfunc->muscle_control=true;
            joint->joint_function->Set_k_p(k_p);
            joint->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(-(T)pi/2,VECTOR<T,3>(1,0,0))));
        }
    }

    bool use_track=parameter_list.Get_Parameter("use_track",false);

    if(num_bodies==5){
        if(parameter_list.Get_Parameter("push_up",false)){
            rigid_body_particles.Rigid_Body(1).frame.r=QUATERNION<T>(3*pi/4,VECTOR<T,3>(0,0,1))*rigid_body_particles.Rigid_Body(1).frame.r;
            rigid_body_particles.Rigid_Body(1).frame.t=VECTOR<T,3>(0,.8,0);

            rigid_body_particles.Rigid_Body(3).Set_Mass(parameter_list.Get_Parameter("mass_scale",(T)100)*rigid_body_particles.Rigid_Body(3).mass);

            arb->joint_mesh.joints(1)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(pi/2,VECTOR<T,3>(1,0,0))));
            arb->joint_mesh.joints(2)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(-3*pi/4,VECTOR<T,3>(1,0,0))));
            arb->joint_mesh.joints(3)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0))));
            arb->joint_mesh.joints(4)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(0,VECTOR<T,3>(1,0,0))));

            //arb->joint_mesh.joints(1)->joint_function->Set_Target_Angle(QUATERNION<T>(pi/2,VECTOR<T,3>(1,0,0)));
            //arb->joint_mesh.joints(2)->joint_function->Set_Target_Angle(QUATERNION<T>(-3*pi/4,VECTOR<T,3>(1,0,0)));
            //arb->joint_mesh.joints(3)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
            //arb->joint_mesh.joints(4)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));

            arb->joint_mesh.joints(1)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));
            arb->joint_mesh.joints(2)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/2,VECTOR<T,3>(1,0,0)));
            arb->joint_mesh.joints(3)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
            arb->joint_mesh.joints(4)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));
        }
        else{
            rigid_body_particles.Rigid_Body(1).frame.r=QUATERNION<T>(3*pi/4,VECTOR<T,3>(0,0,1))*rigid_body_particles.Rigid_Body(1).frame.r;
            rigid_body_particles.Rigid_Body(1).frame.t=VECTOR<T,3>(0,.8,0);

            rigid_body_particles.Rigid_Body(3).Set_Mass(parameter_list.Get_Parameter("mass_scale",(T)100)*rigid_body_particles.Rigid_Body(3).mass);

            arb->joint_mesh.joints(1)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(0,VECTOR<T,3>(1,0,0))));
            arb->joint_mesh.joints(2)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0))));
            arb->joint_mesh.joints(3)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0))));
            arb->joint_mesh.joints(4)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(0,VECTOR<T,3>(1,0,0))));

            arb->joint_mesh.joints(1)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));
            arb->joint_mesh.joints(2)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
            arb->joint_mesh.joints(3)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
            arb->joint_mesh.joints(4)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));

            //arb->joint_mesh.joints(1)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));
            //arb->joint_mesh.joints(2)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/2,VECTOR<T,3>(1,0,0)));
            //arb->joint_mesh.joints(3)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
            //arb->joint_mesh.joints(4)->joint_function->Set_Target_Angle(QUATERNION<T>(0,VECTOR<T,3>(1,0,0)));
        }
    }
    else if(num_bodies==3){
        rigid_body_particles.Rigid_Body(1).frame.r=QUATERNION<T>(pi,VECTOR<T,3>(0,0,1))*rigid_body_particles.Rigid_Body(1).frame.r;
        rigid_body_particles.Rigid_Body(1).frame.t=VECTOR<T,3>(0,.05,0);

        rigid_body_particles.Rigid_Body(2).Set_Mass(parameter_list.Get_Parameter("mass_scale",(T)100)*rigid_body_particles.Rigid_Body(3).mass);

        arb->joint_mesh.joints(1)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(0,VECTOR<T,3>(1,0,0))));
        arb->joint_mesh.joints(2)->Set_Joint_Frame(FRAME_3D<T>(QUATERNION<T>(0,VECTOR<T,3>(1,0,0))));

        arb->joint_mesh.joints(1)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
        arb->joint_mesh.joints(2)->joint_function->Set_Target_Angle(QUATERNION<T>(-pi/4,VECTOR<T,3>(1,0,0)));
    }

    arb->Update_With_Breadth_First_Directed_Graph(1);

    if(parameter_list.Get_Parameter("add_obstacle",false)){
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box",0.3);
        RIGID_BODY<TV>* rigid_body=arb->rigid_body_list.rigid_bodies(id);
        rigid_body->frame.t=VECTOR<T,3>(1,3,0);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("obstacle");
        rigid_body->Set_Mass(1);
        rigid_body->is_static=true;
    }
}
//#####################################################################
// Function Arm
//#####################################################################
void Arm()
{
    //PRISMATIC_JOINT_2D<T>* joint=new PRISMATIC_JOINT_2D<T>();arb->joint_mesh.Add_Joint(joint);
    //joint->Set_Prismatic_Component_Bounds(VECTOR_2D<T>(0,-.5),VECTOR_2D<T>(0,.5));
    std::cout<<" building ARM\n";
    RIGID_BODY<TV> *rigid_body1 = 0;
    RIGID_BODY<TV> *rigid_body2 = 0;

    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box",1);
    rigid_body1=arb->rigid_body_list.rigid_bodies(id);
    rigid_body1->frame.t=VECTOR<T,3>(0,4,0);
    rigid_body1->Set_Coefficient_Of_Restitution(0.5);
    rigid_body1->Set_Coefficient_Of_Friction(0);
    rigid_body1->Set_Name("square1");
    rigid_body1->Set_Mass(1);
    //joint->Set_Joint_To_Parent_Frame(FRAME_3D<T>(VECTOR<T,3>(0,2)));
    
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/box",1);
    rigid_body2=arb->rigid_body_list.rigid_bodies(id);
    rigid_body2->frame.t=VECTOR<T,3>(0,(T)6.7,0);
    rigid_body2->Set_Coefficient_Of_Friction(0);
    rigid_body2->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body2->Set_Name("square2");
    rigid_body2->is_static=true;
    //joint->Set_Joint_To_Child_Frame(FRAME_3D<T>(VECTOR<T,3>(0,-1)));

    muscle=new MUSCLE<T,GRID<TV> >(arb->muscle_list->muscle_force_curve);
    muscle->Set_Attachment_Point_1(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body1,VECTOR<T,3>(0,1,0))); //don't think these are exactly right . . .
    muscle->Set_Attachment_Point_2(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body2,VECTOR<T,3>(0,-1,0))); //don't think these are exactly right . . .
    muscle->Set_Optimal_Length((T).5);
    muscle->Set_Peak_Force(peak_force);
    muscle->Set_Max_Shortening_Velocity(1);
    arb->muscle_list->Add_Muscle(muscle);
/*

    JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(joint,rigid_body1,rigid_body2);
    joint->Set_Joint_Function(jfunc);
    jfunc->Set_Desired_Angle(-.3);
    joint->joint_function->Set_k_p(10);
*/
    //arb->joint_mesh.Add_Articulation(1,2,1);
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    std::cout<<"PREPROCESSING\n"<<std::endl;
    //T activation=fabs(muscle->goal_length/muscle->optimal_length-muscle->Length()/muscle->optimal_length);
    //if(activation>1) activation=1;
    //if(muscle) arb->muscle_list->Set_Muscle_Activation(muscle->id,muscle->Calculate_Activation(muscle2->Force(1)));
    
}
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Write_Output_Files(frame);
#if 0
    const RIGID_BODY_LIST<T,VECTOR<T,3> >& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    VECTOR<T,3> total_linear_momentum,total_angular_momentum;
    for(int i=1;i<=rigid_body_list.Number_Of_Elements();i++){
        LOG::cout << i << ": " << rigid_body_particles.Rigid_Body(i).Angular_Momentum() << std::endl;
        total_linear_momentum+=rigid_body_particles.Rigid_Body(i).mass*rigid_body_particles.Rigid_Body(i).velocity;
        total_angular_momentum+=VECTOR<T,3>::Cross_Product(rigid_body_particles.Rigid_Body(i).frame.t,rigid_body_particles.Rigid_Body(i).mass*rigid_body_particles.Rigid_Body(i).velocity)+rigid_body_particles.Rigid_Body(i).Angular_Momentum();}
    LOG::cout << "MOMENTA === linear " << total_linear_momentum << ", angular " << total_angular_momentum << std::endl;
#endif
}
//#####################################################################
};
}
#endif
