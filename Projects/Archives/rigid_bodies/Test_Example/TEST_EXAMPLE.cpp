//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Guendelman - May 2, 2003
//#####################################################################
#include "TEST_EXAMPLE.h"

using namespace PhysBAM;

template<class T>
TEST_EXAMPLE<T>::TEST_EXAMPLE(int parameter) : RIGID_BODIES_EXAMPLE<T>(parameter)
{
    last_frame=720;
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
    switch (parameter)
    {
        case 0: Dribble_Test(); break;
        case 1: Chute_Test(); break;
    }
}

template<class T> void TEST_EXAMPLE<T>::
Dribble_Test()
{
    VECTOR_3D<T> ground_color=RGB_COLORS<T>().Ground_Tan(),color1=RGB_COLORS<T>().Light_Steel_Blue();

    RIGID_BODY<TV> *rigid_body = 0;

    rigid_body=Initialize_Rigid_Body("subdivided_box");
    rigid_body->Set_Mass(0.1);
    rigid_body->position=VECTOR_3D<T>(0,5,0);
    rigid_body->velocity=VECTOR_3D<T>(0,-2,0);
    rigid_body->Set_Coefficient_Of_Restitution(1.0);
    rigid_body->Set_Name("box");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body("sphere");
    rigid_body->position=VECTOR_3D<T>(0,1,0);
    rigid_body->Set_Coefficient_Of_Restitution(1.0);
    rigid_body->Set_Name("sphere (coeff 1)");
    Append_Diffuse_Hints(color1);

    // PLANE -- no gravity
    rigid_body=Initialize_Rigid_Body("ground");
    rigid_body->Set_Coefficient_Of_Restitution(1.0);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    Append_Diffuse_Hints(ground_color,false);
}

template<class T> void TEST_EXAMPLE<T>::
Chute_Test()
{
    VECTOR_3D<double> ground_color=RGB_COLORS<T>().Ground_Tan(),color1=RGB_COLORS<T>().Light_Steel_Blue();

    RIGID_BODY<TV> *rigid_body = 0;

    T epsilon = 0.0;
    T mu = 100.0;
    T angle1 = 5.625/2.0;
    T angle = pi/2 - angle1*pi/180.0;

#if 1
    rigid_body=Initialize_Rigid_Body("box", 5);
    rigid_body->Set_Mass(5000);
    rigid_body->position=VECTOR_3D<T>(0,150,0);
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(mu);
    rigid_body->Set_Name("box big");
    Append_Diffuse_Hints(color1);
#endif

#if 0
    rigid_body=Initialize_Rigid_Body("box", 0.5);
    rigid_body->Set_Mass(5);
    rigid_body->position=VECTOR_3D<T>(0,50,0);
    rigid_body->Set_Coefficient_Of_Restitution(epsilon);
    rigid_body->Set_Coefficient_Of_Friction(mu);
    rigid_body->Set_Name("box small");
    Append_Diffuse_Hints(color1);
#endif

    rigid_body=Initialize_Rigid_Body("ground", 5);
    rigid_body->orientation=QUATERNION<T>(angle, VECTOR_3D<T>(0,0,1));
    rigid_body->Set_Coefficient_Of_Restitution(1.0);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    Append_Diffuse_Hints(ground_color,false);

    rigid_body=Initialize_Rigid_Body("ground", 5);
    rigid_body->orientation=QUATERNION<T>(-angle, VECTOR_3D<T>(0,0,1));
    rigid_body->Set_Coefficient_Of_Restitution(1.0);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    Append_Diffuse_Hints(ground_color,false);
}

namespace PhysBAM
{
    template TEST_EXAMPLE<double>;
    template TEST_EXAMPLE<float>;
}
