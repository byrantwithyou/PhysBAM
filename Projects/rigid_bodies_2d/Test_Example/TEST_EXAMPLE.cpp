//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "TEST_EXAMPLE.h"

using namespace PhysBAM;

template<class T>
TEST_EXAMPLE<T>::TEST_EXAMPLE(int parameter) : RIGID_BODIES_2D_EXAMPLE<T>(parameter)
{
    last_frame=48;
    output_directory="Test_Example/output";
}

template<class T>
TEST_EXAMPLE<T>::~TEST_EXAMPLE()
{}

//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void TEST_EXAMPLE<T>::
Initialize_Rigid_Bodies()
{
    VECTOR_3D<T> color1=RGB_COLORS<T>().Light_Steel_Blue();

    RIGID_BODY<TV> *rigid_body = 0;

#if 1
    RECTANGULAR_RANDOM_PLACEMENT<T> random_placement(BOX_2D<T>(-20,20,1,100));
    Random_Scene_Generator("square_refined", 100, color1, 1234, random_placement);
    for (int i = 1; i <= rigid_body_parameters.list.rigid_bodies.m; i++)
    {
        rigid_body_parameters.list.rigid_bodies(i)->Set_Coefficient_Of_Restitution(0.5);
    }
#endif
#if 0
    rigid_body=Initialize_Rigid_Body("square");
    rigid_body->position=VECTOR_2D<T>(0,0);
    rigid_body->velocity=VECTOR_2D<T>(0,0);
    rigid_body->orientation=pi/5.0;
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Name("square");
    Append_Diffuse_Hints(color1);
#endif
#if 0 
    rigid_body=Initialize_Rigid_Body("circle",2);
    rigid_body->position=VECTOR_2D<T>(0,120);
    rigid_body->velocity=VECTOR_2D<T>(0,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Name("circle");
    Append_Diffuse_Hints(color1);
#endif

    rigid_body=Initialize_Rigid_Body("ground");
    rigid_body->position=VECTOR_2D<T>(0,-20);
    rigid_body->velocity=VECTOR_2D<T>(0,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    Append_Diffuse_Hints(color1);
}

namespace PhysBAM
{
    template TEST_EXAMPLE<double>;
    template TEST_EXAMPLE<float>;
}
