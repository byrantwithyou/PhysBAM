//#####################################################################
// Copyright 2004-2007, Craig Schroeder, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "ARB_EXAMPLE_3D.h"

using namespace PhysBAM;

template<class T>
ARB_EXAMPLE_3D<T>::ARB_EXAMPLE_3D(const STREAM_TYPE stream_type,int parameter)
    :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solids_parameters)
{
    arb=new ARTICULATED_RIGID_BODY<TV>(solid_body_collection.deformable_object.particles,solids_parameters.rigid_body_parameters.list);
    solids_parameters.rigid_body_parameters.Set_Articulated_Rigid_Body(arb);
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.cfl=(T).1;
    fluids_parameters.simulate=false;
    arb->Set_Iterative_Tolerance((T)1e-4);
    arb->Set_Contact_Level_Iterations(5);
    arb->Set_Shock_Propagation_Level_Iterations(5);
    arb->Set_Poststabilization_Iterations(5);
    arb->Set_Use_Shock_Propagation(false);
    arb->Set_Do_Final_Pass(false);
    last_frame=500;
    frame_rate=24*2;
    output_directory="ARB/output";
    LOG::cout << "Frame rate: "<<frame_rate<<std::endl;
    arb->Set_Iterative_Tolerance((T)1e-6);
    //arb->Set_Constraint_Iterations(100);
    //arb->Set_Collisions_ARB_Iterations(100);
}

template<class T>
ARB_EXAMPLE_3D<T>::~ARB_EXAMPLE_3D()
{}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void ARB_EXAMPLE_3D<T>::
Initialize_Rigid_Bodies()
{
    Coney_Island();
//    Two_Twins_Tumbling();
//    Five_Phallic_Falling();
//    One_Wonder_Wobbling();
//    Diez_Draping_Diddlehoppers();
//    Twelve_Twirps_Twiddling();
}
//#####################################################################
// Initialization Functions
//#####################################################################
template<class T> void ARB_EXAMPLE_3D<T>::
Coney_Island()
{
    JOINT<TV>* joint=new POINT_JOINT<TV>();
    joint->Set_Prismatic_Component_Translation(TV((T)2.5,0,0));
    joint->Use_Twist_Constraint(0,0);
    RIGID_BODY<TV> *box1=0,*box2=0;

    box1=&tests.Add_Rigid_Body("subdivided_box",1,(T).3);
    box1->Set_Coefficient_Of_Restitution((T).5);
    box1->Set_Coefficient_Of_Friction((T)0.5);
    box1->Set_Name("square1");
    box1->is_static=true;
    //box1->Set_Mass(10);
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,1,0),QUATERNION<T>((T)pi/2,TV(0,0,1))));
    //joint2->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(2,1,1)));
    
    box2=&tests.Add_Rigid_Body("subdivided_box",1,(T).3);
    box2->Frame().t=TV(0,(T)4.5,0);
    box2->Angular_Momentum()=TV(4,10,0);//make 0,0,10 with .5 prismatic if you want the stacked example
    // make (10,-3,15) to replicate hinge edge intersection problem with rotation around z
    box2->Set_Coefficient_Of_Restitution((T).5);
    box2->Set_Coefficient_Of_Friction((T)0.5);
    box2->Set_Name("square2");
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-1,0),QUATERNION<T>((T)pi/2,TV(0,0,1))));
    //joint2->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,1,1)));

    arb->joint_mesh.Add_Articulation(1,2,joint);
    //arb->joint_mesh.Add_Articulation(1,2,2);

    //arb->Update_With_Breadth_First_Directed_Graph(1);

    tests.Add_Ground((T).5,-1);
}

template<class T> void ARB_EXAMPLE_3D<T>::
Five_Phallic_Falling()
{
    RIGID_BODY<TV> *rigid_body[5]={0};
    JOINT<TV>* joint[4]={0};
    T x[5]={1,3,5,7,7},y[5]={2,2,2,0,4};

    for(int i=0;i<5;i++){
        rigid_body[i]=&tests.Add_Rigid_Body("subdivided_box",1,(T).3);
        rigid_body[i]->Frame().t=TV(x[i],y[i],0);
        rigid_body[i]->Twist().linear=TV(-5,-14,0);
        rigid_body[i]->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body[i]->Set_Name(STRING_UTILITIES::string_sprintf("square%d",(i+1)));

        if(i>0){joint[i-1]=new POINT_JOINT<TV>;arb->joint_mesh.Add_Articulation(rigid_body[i-1]->particle_index,rigid_body[i]->particle_index,joint[i-1]);
            joint[i-1]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));
            joint[i-1]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,1,1)));}}

    rigid_body[0]->Twist().linear=TV(-5,-14,0);
    rigid_body[3]->Twist().linear=TV(-15,-15,0);
    joint[2]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,-1,1)));
    joint[3]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,1)));

    tests.Add_Ground((T).5,-6);
}
template<class T> void ARB_EXAMPLE_3D<T>::
Two_Twins_Tumbling()
{
    RIGID_BODY<TV> *box1=&tests.Add_Rigid_Body("subdivided_box",1,(T).3);
    box1->Frame().t=TV(0,2,0);
    //box1->Twist().linear=TV(-5,-4,0);
    //box1->Angular_Momentum()=TV(15,0,10);
    box1->Set_Coefficient_Of_Restitution((T)0.5);
    box1->Set_Coefficient_Of_Friction((T)0.5);
    box1->Set_Name("square1");
    
    RIGID_BODY<TV> *box2=&tests.Add_Rigid_Body("subdivided_box",1,(T).3);
    box2->Frame().t=TV(2,2,0);
    box2->Twist().linear=TV(0,0,0);
    box2->Angular_Momentum()=TV(10,-3,20);//make 0,0,10 with .5 prismatic if you want the stacked example
    // make (10,-3,15) to replicate hinge edge intersection problem with rotation around z
    box2->Set_Coefficient_Of_Restitution((T)0.5);
    box2->Set_Coefficient_Of_Friction((T).5);
    box2->Set_Name("square2");

    POINT_JOINT<TV>* joint=new POINT_JOINT<TV>;arb->joint_mesh.Add_Articulation(box1->particle_index,box2->particle_index,joint);
    joint->Use_Twist_Constraint(0,0);
    joint->Use_Theta_Constraint(0,0);

    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));
    //joint2->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(2,1,1)));
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,1,1)));
    //joint2->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,1,1)));

    tests.Add_Ground((T).5,-3);
}
template<class T> void ARB_EXAMPLE_3D<T>::
Three_Things_Throwing()
{

}
template<class T> void ARB_EXAMPLE_3D<T>::
Four_Fools_Flailing()
{

}
template<class T> void ARB_EXAMPLE_3D<T>::
One_Wonder_Wobbling()
{
    RIGID_BODY<TV> *plank=0,*box=0;

    plank=&tests.Add_Rigid_Body("plank",1,(T).3);
    plank->Frame().r=QUATERNION<T>((T)pi/2,TV(1,0,0));
    plank->Set_Coefficient_Of_Restitution((T)0.5);
    plank->Set_Name("static");
    plank->is_static=true;
    
    box=&tests.Add_Rigid_Body("subdivided_box",1,(T).3);
    box->Frame().t=TV(1,6,-(T)1.5);
    box->Frame().r=QUATERNION<T>((T)pi/2,TV(0,0,1));
    box->Twist().linear=TV(-10,0,0);
    box->Set_Coefficient_Of_Restitution((T)0.5);
    box->Set_Name("square2");

    JOINT<TV>* joint=new POINT_JOINT<TV>;
    arb->joint_mesh.Add_Articulation(plank->particle_index,box->particle_index,joint);
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,-(T).5,-5)));
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,1,1)));

    tests.Add_Ground((T).1,-5);
}

template<class T> void ARB_EXAMPLE_3D<T>::
Diez_Draping_Diddlehoppers()
{
    JOINT<TV>* joint[11];
    JOINT<TV>* jointB[11];
    for(int i=0;i<11;i++){joint[i]=new POINT_JOINT<TV>;jointB[i]=new POINT_JOINT<TV>;}

    RIGID_BODY<TV> *rigid_body = 0;
    T scale_factor=(T).65;
    T shift=sqrt((T)2.)/2;

    //green one
    rigid_body=&tests.Add_Rigid_Body("plank",1,(T).3);
    rigid_body->Frame().t=TV(0,3-8*shift,-2-9*shift+(T).25);
    rigid_body->Frame().r=QUATERNION<T>((T)pi/2,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->Set_Name("static");
    rigid_body->is_static=true;
    joint[0]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,(T).5,-5)));
    jointB[0]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,(T).5,-5)));

    //purple one
    rigid_body=&tests.Add_Rigid_Body("plank",1,(T).3);
    rigid_body->Frame().t=TV(0,3-8*shift,2+9*shift-(T).25);
    rigid_body->Frame().r=QUATERNION<T>((T)pi/2,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->Set_Name("static");
    rigid_body->is_static=true;
    joint[10]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-(T).5,-5)));
    jointB[10]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1 ,-(T).5,-5)));
    
    // rope
    // 3
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(0,8-7*shift,-2-7*shift);
    rigid_body->Frame().r=QUATERNION<T>((T)pi/4,TV(1,0,0));
    //rigid_body->Set_Mass(10);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint[0]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,0)));
    joint[1]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,0)));
    jointB[0]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,0)));
    jointB[1]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,0)));
    
    // 4
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(0,8-5*shift,-2-5*shift);
    rigid_body->Frame().r=QUATERNION<T>((T)pi/4,TV(1,0,0));
    //rigid_body->Set_Mass(10);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint[1]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,0)));
    joint[2]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,0)));
    jointB[1]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,0)));
    jointB[2]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,0)));

    // 5
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(0,8-3*shift,-2-3*shift);
    rigid_body->Frame().r=QUATERNION<T>((T)pi/4,TV(1,0,0));
    //rigid_body->Set_Mass(10);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint[2]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,0)));
    joint[3]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,0)));
    jointB[2]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,0)));
    jointB[3]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,0)));

    // 6
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(0,8-shift,-2-shift);
    rigid_body->Frame().r=QUATERNION<T>((T)pi/4,TV(1,0,0));
    //rigid_body->Set_Mass(10);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint[3]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,0)));
    joint[4]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,0)));
    jointB[3]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,0)));
    jointB[4]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,0)));
    
    // 7
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(0,8,-1);
    rigid_body->Frame().r=QUATERNION<T>((T)pi/2,TV(1,0,0));
    //rigid_body->Set_Mass(10);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint[4]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,0)));
    joint[5]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,0)));
    jointB[4]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,0)));
    jointB[5]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,0)));

    // 8
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(0,8,1);
    rigid_body->Frame().r=QUATERNION<T>((T)pi/2,TV(1,0,0));
    rigid_body->Set_Mass(10);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint[5]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,0)));
    joint[6]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,0)));
    jointB[5]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,0)));
    jointB[6]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,0)));


    // 9
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(0,8-shift,2+shift);
    rigid_body->Frame().r=QUATERNION<T>(3*(T)pi/4,TV(1,0,0));
    rigid_body->Set_Mass(10);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint[6]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,0)));
    joint[7]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,0)));
    jointB[6]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,0)));
    jointB[7]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,0)));

    // 10
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(0,8-3*shift,2+3*shift);
    rigid_body->Frame().r=QUATERNION<T>(3*(T)pi/4,TV(1,0,0));
    //rigid_body->Set_Mass(10);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint[7]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,0)));
    joint[8]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,0)));
    jointB[7]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,0)));
    jointB[8]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,0)));

    // 11
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(0,8-5*shift,2+5*shift);
    rigid_body->Frame().r=QUATERNION<T>(3*(T)pi/4,TV(1,0,0));
    //rigid_body->Set_Mass(10);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint[8]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,0)));
    joint[9]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,0)));
   jointB[8]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,0)));
    jointB[9]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,0)));

    // 12
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(0,8-7*shift,2+7*shift);
    rigid_body->Frame().r=QUATERNION<T>(3*(T)pi/4,TV(1,0,0));
    //rigid_body->Set_Mass(10);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint[9]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,0)));
    joint[10]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,0)));
    jointB[9]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,0)));
    jointB[10]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,0)));

    arb->joint_mesh.Add_Articulation(1,3,1);
    arb->joint_mesh.Add_Articulation(1,3,2);
    arb->joint_mesh.Add_Articulation(3,4,3);
    arb->joint_mesh.Add_Articulation(3,4,4);
    arb->joint_mesh.Add_Articulation(4,5,5);
    arb->joint_mesh.Add_Articulation(4,5,6);
    arb->joint_mesh.Add_Articulation(5,6,7);
    arb->joint_mesh.Add_Articulation(5,6,8);
    arb->joint_mesh.Add_Articulation(6,7,9);
    arb->joint_mesh.Add_Articulation(6,7,10);
    arb->joint_mesh.Add_Articulation(7,8,11);
    arb->joint_mesh.Add_Articulation(7,8,12);
    arb->joint_mesh.Add_Articulation(8,9,13);
    arb->joint_mesh.Add_Articulation(8,9,14);
    arb->joint_mesh.Add_Articulation(9,10,15);
    arb->joint_mesh.Add_Articulation(9,10,16);
    arb->joint_mesh.Add_Articulation(10,11,17);
    arb->joint_mesh.Add_Articulation(10,11,18);
    arb->joint_mesh.Add_Articulation(11,12,19);
    arb->joint_mesh.Add_Articulation(11,12,20);
    arb->joint_mesh.Add_Articulation(12,2,21);
    arb->joint_mesh.Add_Articulation(12,2,22);
    /*
    arb->joint_mesh.Add_Articulation(1,3,1);
    arb->joint_mesh.Add_Articulation(3,4,2);
    arb->joint_mesh.Add_Articulation(4,5,3);
    arb->joint_mesh.Add_Articulation(5,6,4);
    arb->joint_mesh.Add_Articulation(6,7,5);
    arb->joint_mesh.Add_Articulation(7,8,6);
    arb->joint_mesh.Add_Articulation(8,9,7);
    arb->joint_mesh.Add_Articulation(9,10,8);
    arb->joint_mesh.Add_Articulation(10,11,9);
    arb->joint_mesh.Add_Articulation(11,12,10);
    arb->joint_mesh.Add_Articulation(12,2,11);
*/
    tests.Add_Ground((T).12,-5);
}

template<class T> void ARB_EXAMPLE_3D<T>::
Twelve_Twirps_Twiddling()
{
    JOINT<TV>* joint=arb->joint_mesh.joints(arb->joint_mesh.Add_Joint());
    JOINT<TV>* joint2=arb->joint_mesh.joints(arb->joint_mesh.Add_Joint());
    JOINT<TV>* joint3=arb->joint_mesh.joints(arb->joint_mesh.Add_Joint());
    JOINT<TV>* joint4=arb->joint_mesh.joints(arb->joint_mesh.Add_Joint());
    JOINT<TV>* joint5=arb->joint_mesh.joints(arb->joint_mesh.Add_Joint());
    JOINT<TV>* joint6=arb->joint_mesh.joints(arb->joint_mesh.Add_Joint());
    JOINT<TV>* joint7=arb->joint_mesh.joints(arb->joint_mesh.Add_Joint());
    JOINT<TV>* joint8=arb->joint_mesh.joints(arb->joint_mesh.Add_Joint());
    JOINT<TV>* joint9=arb->joint_mesh.joints(arb->joint_mesh.Add_Joint());
    JOINT<TV>* joint10=arb->joint_mesh.joints(arb->joint_mesh.Add_Joint());
    JOINT<TV>* joint11=arb->joint_mesh.joints(arb->joint_mesh.Add_Joint());
    JOINT<TV>* joint12=arb->joint_mesh.joints(arb->joint_mesh.Add_Joint());

    RIGID_BODY<TV> *rigid_body = 0;
    T scale_factor=1;
    T shift=sqrt((T)2.)/2;

    // 1
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(0,4,0);
    rigid_body->Twist().linear=TV(-10,0,0);
    //rigid_body->Frame().r=QUATERNION<T>((T)pi/4,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,1)));
    joint2->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,-1,1)));
    
    // 2
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(2,4,0);
    rigid_body->Twist().linear=TV(-10,0,0);
    //rigid_body->Frame().r=QUATERNION<T>((T)pi/4,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint2->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,1)));
    joint3->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,-1,1)));

    // 3
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(4,2,0);
    rigid_body->Twist().linear=TV(0,10,0);
    //rigid_body->Frame().r=QUATERNION<T>((T)pi/4,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint3->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,1,1)));
    joint4->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,-1,1)));

    // 4
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(4,0,0);
    rigid_body->Twist().linear=TV(0,10,0);
    //rigid_body->Frame().r=QUATERNION<T>((T)pi/4,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint4->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,1,1)));
    joint5->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,-1,1)));
    
    // 5
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(4,-2,0);
    rigid_body->Twist().linear=TV(0,10,0);
//    rigid_body->Frame().r=QUATERNION<T>((T)pi/2,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint5->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,1,1)));
    joint6->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,-1,1)));

    // 6
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(2,-4,0);
    rigid_body->Twist().linear=TV(10,0,0);
    //rigid_body->Frame().r=QUATERNION<T>((T)pi/2,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint6->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,1,1)));
    joint7->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,1)));

    // 7
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(0,-4,0);
    rigid_body->Twist().linear=TV(10,0,0);
    //rigid_body->Frame().r=QUATERNION<T>(3*(T)pi/4,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint7->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,1,1)));
    joint8->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,1)));

    // 8
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(-2,-4,0);
    rigid_body->Twist().linear=TV(10,0,0);
    //rigid_body->Frame().r=QUATERNION<T>(3*(T)pi/4,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint8->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,1,1)));
    joint9->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1,1)));

    // 9
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(-4,-2,0);
    rigid_body->Twist().linear=TV(0,-10,0);
    //rigid_body->Frame().r=QUATERNION<T>(3*(T)pi/4,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint9->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,1)));
    joint10->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));

    // 10
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(-4,0,0);
    rigid_body->Twist().linear=TV(0,-10,0);
    //rigid_body->Frame().r=QUATERNION<T>(3*(T)pi/4,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint10->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,1)));
    joint11->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));

    // 11
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(-4,2,0);
    rigid_body->Twist().linear=TV(0,-10,0);
    //rigid_body->Frame().r=QUATERNION<T>(3*(T)pi/4,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint11->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,1)));
    joint12->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));

    // 12
    rigid_body=&tests.Add_Rigid_Body("subdivided_box",scale_factor,(T).3);
    rigid_body->Frame().t=TV(-2,4,0);
    rigid_body->Twist().linear=TV(-10,0,20);
    //rigid_body->Frame().r=QUATERNION<T>(3*(T)pi/4,TV(1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    joint12->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1,1)));
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,-1,1)));

    arb->joint_mesh.Add_Articulation(12,1,1);
    arb->joint_mesh.Add_Articulation(1,2,2);
    arb->joint_mesh.Add_Articulation(2,3,3);
    arb->joint_mesh.Add_Articulation(3,4,4);
    arb->joint_mesh.Add_Articulation(4,5,5);
    arb->joint_mesh.Add_Articulation(5,6,6);
    arb->joint_mesh.Add_Articulation(6,7,7);
    arb->joint_mesh.Add_Articulation(7,8,8);
    arb->joint_mesh.Add_Articulation(8,9,9);
    arb->joint_mesh.Add_Articulation(9,10,10);
    arb->joint_mesh.Add_Articulation(10,11,11);
    arb->joint_mesh.Add_Articulation(11,12,12);

    tests.Add_Ground((T).12,-20);
}

namespace PhysBAM
{
    template class ARB_EXAMPLE_3D<double>;
    template class ARB_EXAMPLE_3D<float>;
}
