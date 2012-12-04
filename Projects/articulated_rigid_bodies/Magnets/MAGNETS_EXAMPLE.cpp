//#####################################################################
// Copyright 2004-2007, Craig Schroeder, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/PRISMATIC_TWIST_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include "MAGNETS_EXAMPLE.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> MAGNETS_EXAMPLE<T>::
MAGNETS_EXAMPLE(const STREAM_TYPE stream_type)
    :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection)
{
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    last_frame=500;
    frame_rate=24*4;
    std::cout<<"Frame rate: "<<frame_rate<<std::endl;

    write_last_frame=true;

    rg.Set_Seed(2118);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> MAGNETS_EXAMPLE<T>::
~MAGNETS_EXAMPLE()
{
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class T> void MAGNETS_EXAMPLE<T>::
Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class T> void MAGNETS_EXAMPLE<T>::
Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    output_directory="Magnets/output";
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T> void MAGNETS_EXAMPLE<T>::
Initialize_Bodies()
{
    RIGID_BODY<TV> *rigid_body=0;

    arb=&solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb->Set_Iterative_Tolerance((T)1e-4);
//    arb->Set_Extra_Iterations_Per_Contact_Level_Factor(100);
//    arb->Set_Extra_Iterations_Per_Shock_Propagation_Level_Factor(100);
    arb->Set_Poststabilization_Iterations(100);
    //arb->Set_Use_Shock_Propagation(false);
    arb->Set_Do_Final_Pass(false);

    T up_shift=31;
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",11,(T).5);
    rigid_body->Frame().t=TV(0,0,(T)11.5);
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->name="parent";
    rigid_body->is_static=true;
    parent_id=rigid_body->particle_index;
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",11,(T).5);
    rigid_body->Frame().t=TV(0,22,(T)11.5);
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->name="parent";
    rigid_body->is_static=true;
    
    tests.Add_Ground((T).5,-11,0);

    Random_Letters(20,-9); // 20 for full fridge
    Drop_Letter("S_xsub",parent_id,TV(-7,up_shift,0),ROTATION<TV>(),true);//pi/8
    Drop_Letter("I_xsub",parent_id,TV(-5,up_shift,0),ROTATION<TV>((T)(T)pi/12,TV(0,0,1)),true);
    Drop_Letter("G_xsub",parent_id,TV(-3,up_shift,0),ROTATION<TV>((T)(T)pi/6,TV(0,0,1)),true);//pi/10
    Drop_Letter("G_xsub",parent_id,TV(-1,up_shift,0),ROTATION<TV>(),true);//pi/10
    Drop_Letter("R_xsub",parent_id,TV(1,up_shift,0),ROTATION<TV>(),true);//pi/4
    Drop_Letter("A_xsub",parent_id,TV(3,up_shift,0),ROTATION<TV>((T)(T)pi/12,TV(0,0,1)),true);//pi/64
    Drop_Letter("P_xsub",parent_id,TV(5,up_shift,0),ROTATION<TV>(),true);//pi/10
    Drop_Letter("H_xsub",parent_id,TV(7,up_shift,0),ROTATION<TV>((T)(T)pi/15,TV(0,0,1)),true);
    Drop_Letter("A_xsub",parent_id,TV(-9,up_shift,0),ROTATION<TV>(-(T)(T)pi/3,TV(0,0,1)),false);
    Drop_Letter("G_xsub",parent_id,TV(9,up_shift,0),ROTATION<TV>(),false);

    rigid_body=&tests.Add_Rigid_Body("plank",(T).5,(T).5);
    rigid_body->Frame().t=TV(10,17+(T)2.375,(T)1.5);
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->name="handle";
    rigid_body->is_static=true;
    rigid_body=&tests.Add_Rigid_Body("plank",(T).5,(T).5);
    rigid_body->Frame().t=TV(10,17-(T)2.375,(T)1.5);
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->name="handle";
    rigid_body->is_static=true;
    rigid_body=&tests.Add_Rigid_Body("plank",(T).5,(T).5);
    rigid_body->Frame().t=TV(10,17,-1);
    rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->name="handle";
    rigid_body->is_static=true;
    
    // add forces
    tests.Add_Gravity();
    solid_body_collection.Update_Simulated_Particles();

    BASE::Initialize_Bodies();
}
//#####################################################################
// Function Random_Letters
//#####################################################################
template<class T> void MAGNETS_EXAMPLE<T>::
Random_Letters(int levels,T up_shift)
{
    for(int i=0;i<levels;i++) for(int j=-4;j<6;j++){
        T x_shift=rg.Get_Uniform_Number(-(T).24,(T).24),y_shift=rg.Get_Uniform_Number(-(T).24,(T).24),o_shift=rg.Get_Uniform_Number(-(T)3.14,(T)3.14);
        std::string letter;
        if(((i==11 || i==12)&&j==5) || (i==7&&j==3) || (i==11&&j==2)) letter="I_xsub"; //TODO: why was this called here? rg.Get_Uniform_Integer(1,7);
        else if((i==3 && j==2) || (i==5&&j==-1)) letter="R_xsub"; // rg.Get_Uniform_Integer(1,7);
        else if((j==-4 && (i==3 || i==6)) || (i==8&&j==1) || (i==7&&j==-1)) letter="H_xsub"; // rg.Get_Uniform_Integer(1,7);
        else letter=Get_Random_Letter();
        Drop_Letter(letter,parent_id,TV(j*2-1+x_shift,up_shift+2*i+y_shift,0),ROTATION<TV>(o_shift,TV(0,0,1)),false);}
}
//#####################################################################
// Function Get_Random_Letter
//#####################################################################
template<class T> std::string MAGNETS_EXAMPLE<T>::
Get_Random_Letter()
{
    int letter=rg.Get_Uniform_Integer(1,7);
    switch(letter){
        case 1: return "S_xsub";
        case 2: return "I_xsub";
        case 3: return "G_xsub";
        case 4: return "R_xsub";
        case 5: return "A_xsub";
        case 6: return "P_xsub";
        case 7: default: return "H_xsub";}
}
//#####################################################################
// Function Drop_Letter
//#####################################################################
template<class T> void MAGNETS_EXAMPLE<T>::
Drop_Letter(std::string letter,int parent_id,TV start,ROTATION<TV> orient,bool stop_in_middle)
{
    RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("Letters/"+letter,(T)1.9,(T).5);
    rigid_body.Frame().t=start;
    rigid_body.Frame().r=orient*ROTATION<TV>((T)pi,TV(0,1,0));
    rigid_body.Set_Coefficient_Of_Restitution(0);
    rigid_body.name="letter";

    PRISMATIC_TWIST_JOINT<TV>* joint=new PRISMATIC_TWIST_JOINT<TV>;
    arb->joint_mesh.Add_Articulation(parent_id,rigid_body.particle_index,joint);
    if(stop_in_middle) joint->Set_Prismatic_Constraints(VECTOR<bool,3>(true,true,true),TV(0,-8,-10),TV(0,1000,10));
    else joint->Set_Prismatic_Constraints(VECTOR<bool,3>(true,true,true),TV(0,-1000,-10),TV(0,1000,10));
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,22,-(T)11.5),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0))));
}
//#####################################################################
namespace PhysBAM{
template class MAGNETS_EXAMPLE<float>;
template class MAGNETS_EXAMPLE<double>;
}
