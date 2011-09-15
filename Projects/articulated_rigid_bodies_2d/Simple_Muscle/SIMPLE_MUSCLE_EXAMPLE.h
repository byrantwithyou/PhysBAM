//#####################################################################
// Copyright 2004-2005, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLE_MUSCLE_EXAMPLE
//##################################################################### 
#ifndef __SIMPLE_MUSCLE_EXAMPLE__
#define __SIMPLE_MUSCLE_EXAMPLE__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <fstream>
#include "../ARB_PARAMETERS.h"
#include <Forces_And_Torques/ARB_SPRING_CONSTRAINT.h>
namespace PhysBAM{

template<class T,class RW>
class SIMPLE_MUSCLE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW> BASE;
    typedef CONSTRAINED_POINT_IN_RIGID_BODY<T,VECTOR_2D<T>,RIGID_BODY<TV> > T_CONSTRAINED_POINT_IN_RIGID_BODY;

    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;
    using BASE::restart;using BASE::restart_frame;using BASE::output_directory;
    using BASE::solids_parameters;using BASE::write_last_frame;using BASE::data_directory;

    ARTICULATED_RIGID_BODY<TV>* arb;
    MUSCLE<T,GRID<TV> >* muscle;
    ARB_SPRING_CONSTRAINT<T,GRID<TV> > *muscle2;
    bool add_ground;
    PARAMETER_LIST parameter_list;

    SIMPLE_MUSCLE_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>(0,FLUIDS_PARAMETERS_2D<T>::NONE)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=.1;
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

        arb->Use_Muscle_Actuators();
        if(getenv("USE_PD")) arb->Use_PD_Actuators();

        muscle=0;
        
        ARB_PARAMETERS::Read_Common_Parameters("Simple_Muscle/example.param",*this,parameter_list);
    }

//#####################################################################
// Function Output_Curves -- for debugging
//#####################################################################
void Output_Curve(const GRID<TV>& grid,const ARRAY<T,VECTOR<int,1> >& values,const std::string& filename)
{std::ofstream output(filename.c_str());for(int i=1;i<=grid.m;i++) output<<grid.x(i)<<" "<<values(i)<<std::endl;}
void Output_Curves()
{
    MUSCLE_FORCE_CURVE<T>& muscle_force_curve=arb->muscle_list->muscle_force_curve;
    Output_Curve(muscle_force_curve.passive_force_grid,muscle_force_curve.passive_force,"test_passive_force");
    Output_Curve(muscle_force_curve.passive_force_slope_grid,muscle_force_curve.passive_force_slope,"test_passive_force_slope");
    Output_Curve(muscle_force_curve.tendon_force_grid,muscle_force_curve.tendon_force,"test_tendon_force");
    Output_Curve(muscle_force_curve.tendon_force_slope_grid,muscle_force_curve.tendon_force_slope,"test_tendon_force_slope");
    Output_Curve(muscle_force_curve.tendon_length_grid,muscle_force_curve.tendon_length,"test_tendon_length");
    Output_Curve(muscle_force_curve.active_force_grid,muscle_force_curve.active_force,"test_active_force");
    Output_Curve(muscle_force_curve.active_force_slope_grid,muscle_force_curve.active_force_slope,"test_active_force_slope");
    Output_Curve(muscle_force_curve.velocity_grid,muscle_force_curve.velocity_curve,"test_velocity_curve");
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    arb->muscle_list->muscle_force_curve.Initialize(data_directory); // initialize here rather than constructor since data directory might be set after constructor
    //Output_Curves();

    add_ground=true;
    //Arm();
    //Spring();
    Simple_Muscle_Across_Joint();

    if(add_ground){
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.1);
        RIGID_BODY<TV>* rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
        rigid_body->frame.t=VECTOR_2D<T>(0,0);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("ground");
        rigid_body->is_static=true;
        rigid_body->add_to_spatial_partition=false;}

    std::cout <<"\nnumber of muscles is: "<<arb->muscle_list->muscles.m<<std::endl;

    RIGID_BODY_LIST_2D<T>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    for(int i=1;i<=rigid_body_list.Number_Of_Elements();i++) if(!rigid_body_particles.Rigid_Body(i).is_static)
       rigid_body_particles.Rigid_Body(i).Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,(T)0);

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
// Function Simple_Muscle_Across_Joint
//#####################################################################
void Simple_Muscle_Across_Joint()
{
    add_ground=false;
    solids_parameters.gravity=0;

    RIGID_BODY<TV>* rigid_body=0,*parent=0,*child=0;
    POINT_JOINT<TV>* joint=new POINT_JOINT<TV>();
    MUSCLE<T,GRID<TV> > *flexor=new MUSCLE<T,GRID<TV> >(arb->muscle_list->muscle_force_curve);flexor->name="flexor";
    MUSCLE<T,GRID<TV> > *extensor=new MUSCLE<T,GRID<TV> >(arb->muscle_list->muscle_force_curve);extensor->name="extensor";

    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01,true,false,false);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body->frame.r=COMPLEX<T>::Unit_Polar(pi/2);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0);
    rigid_body->Set_Name("parent");
    rigid_body->Set_Mass(1);
    parent=rigid_body;
    joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(-1,0)));

    flexor->Set_Attachment_Point_1(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body,VECTOR_2D<T>(0.3,0)));

    extensor->Set_Attachment_Point_1(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body,VECTOR_2D<T>(-0.5,0)));
    extensor->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body,VECTOR_2D<T>(-1,0.05)));
    
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground",.01,true,false,false);
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Coefficient_Of_Friction(0);
    rigid_body->Set_Name("child");
    rigid_body->Set_Mass(1);
    child=rigid_body;
    joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(-1,0)));

    flexor->Set_Attachment_Point_2(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body,VECTOR_2D<T>(-0.2,0)));
    //flexor->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body,VECTOR_2D<T>(-0.8,0.05)));

    extensor->Set_Attachment_Point_2(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body,VECTOR_2D<T>(-0.5,0)));
    extensor->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body,VECTOR_2D<T>(-1,-0.05)));

    int joint_id=arb->joint_mesh.Add_Joint(joint);
    arb->Add_Articulation(id-1,id,joint_id);

    JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(joint,parent,child);
    joint->Set_Joint_Function(jfunc);

    jfunc->Set_Target_Angle(-(T)1.5);
#if 0
    INTERPOLATION_CURVE<T,T>* angle=new INTERPOLATION_CURVE<T,T>();
    angle->Add_Control_Point(0,-1);
    angle->Add_Control_Point(5,-1);
    angle->Add_Control_Point(10,-1.5);
    angle->Add_Control_Point(15,-1.5);
    angle->Add_Control_Point(20,-0.5);
    angle->Add_Control_Point(25,-1);
    jfunc->interpolation_curve=angle;
#endif

    joint->joint_function->Set_k_p(10);

    joint->Set_Joint_Frame(FRAME_2D<T>(COMPLEX<T>::Unit_Polar(-(T)1)));
    arb->Update_With_Breadth_First_Directed_Graph(1);

    T total_length=flexor->Total_Length();
    flexor->Set_Optimal_Length((T).8*total_length);
    flexor->Set_Tendon_Slack_Length((T).2*total_length);
    flexor->Set_Peak_Force(10);
    flexor->Set_Max_Shortening_Velocity(1);
    arb->muscle_list->Add_Muscle(flexor);

    total_length=extensor->Total_Length();
    extensor->Set_Optimal_Length((T).8*total_length);
    extensor->Set_Tendon_Slack_Length((T).2*total_length);
    extensor->Set_Peak_Force(10);
    extensor->Set_Max_Shortening_Velocity(1);
    arb->muscle_list->Add_Muscle(extensor);
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

    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/square_refined",1,true,false,false);
    rigid_body1=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body1->frame.t=VECTOR_2D<T>(0,4);
    rigid_body1->Set_Coefficient_Of_Restitution(0.5);
    rigid_body1->Set_Coefficient_Of_Friction(0);
    rigid_body1->Set_Name("square1");
    rigid_body1->Set_Mass(1);
    //joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(0,2)));
    
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/square_refined",1,true,false,false);
    rigid_body2=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body2->frame.t=VECTOR_2D<T>(0,6.7);
    rigid_body2->Set_Coefficient_Of_Restitution(0.5);
    rigid_body2->Set_Coefficient_Of_Friction(0);
    rigid_body2->Set_Name("square2");
    rigid_body2->is_static=true;
    //joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(0,-1)));

    muscle=new MUSCLE<T,GRID<TV> >(arb->muscle_list->muscle_force_curve);
    muscle->Set_Attachment_Points(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body1,VECTOR_2D<T>(0,1)),new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body2,VECTOR_2D<T>(0,-1))); //don't think these are exactly right . . .
    muscle->Set_Optimal_Length(.5);
    muscle->Set_Peak_Force(10);
    muscle->Set_Max_Shortening_Velocity(1);
//    muscle->goal_length=.5;
    arb->muscle_list->Add_Muscle(muscle);
/*

    JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(joint,rigid_body1,rigid_body2);
    joint->Set_Joint_Function(jfunc);
    jfunc->Set_Desired_Angle(-.3);
    joint->joint_function->Set_k_p(10);
*/
    //arb->Add_Articulation(1,2,1);
}
//#####################################################################
// Function Spring
//#####################################################################
void Spring()
{
    //PRISMATIC_JOINT_2D<T>* joint=new PRISMATIC_JOINT_2D<T>();arb->joint_mesh.Add_Joint(joint);
    //joint->Set_Prismatic_Component_Bounds(VECTOR_2D<T>(0,-.5),VECTOR_2D<T>(0,.5));
    std::cout<<" building spring\n";
    RIGID_BODY<TV> *rigid_body1 = 0;
    RIGID_BODY<TV> *rigid_body2 = 0;

    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/square_refined",1,true,false,false);
    rigid_body1=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body1->frame.t=VECTOR_2D<T>(4,4);
    rigid_body1->Set_Coefficient_Of_Restitution(0.5);
    rigid_body1->Set_Coefficient_Of_Friction(0);
    rigid_body1->Set_Name("square1");
    rigid_body1->Set_Mass(1);
    //joint->Set_Joint_To_Parent_Frame(FRAME_2D<T>(VECTOR_2D<T>(0,2)));
    
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/square_refined",1,true,false,false);
    rigid_body2=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body2->frame.t=VECTOR_2D<T>(4,6.7);
    rigid_body2->Set_Coefficient_Of_Restitution(0.5);
    rigid_body2->Set_Coefficient_Of_Friction(0);
    rigid_body2->Set_Name("square2");
    rigid_body2->is_static=true;
    //joint->Set_Joint_To_Child_Frame(FRAME_2D<T>(VECTOR_2D<T>(0,-1)));
    
    muscle2=new ARB_SPRING_CONSTRAINT<T,GRID<TV> >();
    muscle2->Set_Attachment_Points(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body1,VECTOR_2D<T>(0,1)),new T_CONSTRAINED_POINT_IN_RIGID_BODY(*rigid_body2,VECTOR_2D<T>(0,-1))); //don't think these are exactly right . . .
    muscle2->Set_Optimal_Length(.5);
    muscle2->target_length=.5;
    arb->muscle_list->Add_Muscle(muscle2);
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
//#####################################################################
};
}
#endif
