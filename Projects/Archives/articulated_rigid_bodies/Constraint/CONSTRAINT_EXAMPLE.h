//#####################################################################
// Copyright 2004, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSTRAINT_EXAMPLE
//#####################################################################
#ifndef CONSTRAINT_EXAMPLE__
#define CONSTRAINT_EXAMPLE__

#include <iostream>
#include "ARTICULATED_RB_EXAMPLE.h"
#include "GRAPH.h"
#include "JOINT.h"
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NR0.h>

namespace PhysBAM {

template<class T>
class CONSTRAINT_EXAMPLE : public ARTICULATED_RB_EXAMPLE<T> {

public:
    //Variables
    GRAPH<T> node_graph;
    T increment;
    T inc1,inc2,inc3;
    bool test1,test2,test3;
    bool toggle1,toggle2,toggle3;
    RANDOM_NR0 rand_num;
    GRAPH_NODE<T>* gnode;
    GRAPH_NODE<T>* test_child;
    GRAPH_NODE<T>* test_child2;
    GRAPH_NODE<T>* end_effector;

    //Methods
    CONSTRAINT_EXAMPLE() {
    strcpy(example_name,"CONSTRAINT Example");
    strcpy(output_directory,"CONSTRAINT_Example");
    initial_time=0;
    final_time=3;
    frame_rate=24;
    increment = 0.05;
    
    inc1=rand_num.Get_Number();
    inc2=rand_num.Get_Number();
    inc3=rand_num.Get_Number();
    }

    ~CONSTRAINT_EXAMPLE(){
    //TODO
    }

    void Initialize() {
    gnode = new GRAPH_NODE<T>(&node_graph, new JOINT<T>());
    node_graph.Set_Root(gnode);

    RIGID_BODY<TV>* rigid_body;
    rigid_body=Initialize_Rigid_Body("../../Public_Data/Rigid_Bodies/plank",0.2);
    rigid_body->Set_Name("test_cylinder1");
    rigid_body->is_static=true; // not sure what this does yet . . . taking from Eran's code
    rigid_body->position=VECTOR_3D<T>(0,0,0);
    rigid_body->Set_Coefficient_Of_Friction(0.6);    
    ((JOINT<T>*)gnode->node_item)->child=rigid_body;
    ((JOINT<T>*)gnode->node_item)->Set_Child_2_Joint_Frame(FRAME<T>(VECTOR_3D<T>(0,0,-1)));

    std::cout << "okay\n";
    test_child = gnode->Add_Child(new JOINT<T>());
    ((JOINT<T>*)test_child->node_item)->parent=rigid_body;
    rigid_body=Initialize_Rigid_Body("../../Public_Data/Rigid_Bodies/plank",0.2);
    rigid_body->Set_Name("test_cylinder2");
    //rigid_body->position=VECTOR_3D<T>(0,-2,0);
    rigid_body->Set_Coefficient_Of_Friction(0.6);
    ((JOINT<T>*)test_child->node_item)->child=rigid_body;
    ((JOINT<T>*)test_child->node_item)->Set_Joint_2_Parent_Frame(FRAME<T>(VECTOR_3D<T>(0,0,-1)));
    ((JOINT<T>*)test_child->node_item)->Set_Child_2_Joint_Frame(FRAME<T>(VECTOR_3D<T>(0,0,-1)));

    test_child2 = test_child->Add_Child(new JOINT<T>());
    ((JOINT<T>*)test_child2->node_item)->parent = rigid_body;
    rigid_body=Initialize_Rigid_Body("../../Public_Data/Rigid_Bodies/plank",0.2);
    rigid_body->Set_Name("test_cylinder3");
    rigid_body->position=VECTOR_3D<T>(0,2,0);
    rigid_body->Set_Coefficient_Of_Friction(0.6);
    ((JOINT<T>*)test_child2->node_item)->child=rigid_body;
    ((JOINT<T>*)test_child2->node_item)->Set_Joint_2_Parent_Frame(FRAME<T>(VECTOR_3D<T>(0,0,-1)));
    ((JOINT<T>*)test_child2->node_item)->Set_Child_2_Joint_Frame(FRAME<T>(VECTOR_3D<T>(0,0,-1)));

    end_effector = test_child2->Add_Child(new JOINT<T>());
    ((JOINT<T>*)end_effector->node_item)->parent = rigid_body;
    rigid_body=Initialize_Rigid_Body("../../Public_Data/Rigid_Bodies/sphere",0.2);
    rigid_body->Set_Name("end_effector");
    rigid_body->Set_Coefficient_Of_Friction(0.6);
    ((JOINT<T>*)end_effector->node_item)->child = rigid_body;
    ((JOINT<T>*)end_effector->node_item)->Set_Joint_2_Parent_Frame(FRAME<T>(VECTOR_3D<T>(0,0,-1)));
    
    ((JOINT<T>*)gnode->node_item)->Set_Constraints(-1,1,-1,1,-1,1);
    ((JOINT<T>*)test_child->node_item)->Set_Constraints(-1,1,-1,1,-1,1);
    ((JOINT<T>*)test_child2->node_item)->Set_Constraints(-1,1,-1,1,-1,1);
    }

    void Run() {
    

    test1 = ((JOINT<T>*)gnode->node_item)->Set_Joint_Frame(FRAME<T>(QUATERNION<T>(inc1,VECTOR_3D<T>(1,0,0))));
    test2 = ((JOINT<T>*)test_child->node_item)->Set_Joint_Frame(FRAME<T>(QUATERNION<T>(inc2,VECTOR_3D<T>(0,1,0))));
    test3 = ((JOINT<T>*)test_child2->node_item)->Set_Joint_Frame(FRAME<T>(QUATERNION<T>(inc3,VECTOR_3D<T>(0,0,1))));
    node_graph.Update_Forward();

    if (!test1) {if (toggle1) toggle1=false; else toggle1=true;}
    if (!test2) {if (toggle2) toggle2=false; else toggle2=true;}
    if (!test3) {if (toggle3) toggle3=false; else toggle3=true;}

    if (toggle1) inc1+=increment; else inc1-=increment;
    if (toggle2) inc2+=increment; else inc2-=increment;
    if (toggle3) inc3+=increment; else inc3-=increment;
    }
};
}

#endif
