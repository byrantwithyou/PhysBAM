//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Guendelman - January 5, 2003
//#####################################################################
#include "ROLLING_FRICTION_EXAMPLE.h"

using namespace PhysBAM;

template<class T>
ROLLING_FRICTION_EXAMPLE<T>::ROLLING_FRICTION_EXAMPLE(int parameter) : RIGID_BODIES_EXAMPLE<T>(parameter)
{
    last_frame=720;
    output_directory="Rolling_Friction_Example/output";
}

template<class T>
ROLLING_FRICTION_EXAMPLE<T>::~ROLLING_FRICTION_EXAMPLE()
{}

//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void ROLLING_FRICTION_EXAMPLE<T>::
Initialize_Rigid_Bodies()
{
    VECTOR_3D<double> ground_color=RGB_COLORS<T>().Ground_Tan(),color1=RGB_COLORS<T>().Light_Steel_Blue();

    RIGID_BODY<TV> *rigid_body = 0;

    if (parameter == 0)
    {
        rigid_body=Initialize_Rigid_Body("sphere");
        rigid_body->position=VECTOR_3D<T>(-2,1,45);
        rigid_body->velocity=VECTOR_3D<T>(0,0,-10);
        rigid_body->Set_Coefficient_Of_Friction(0.1);
        rigid_body->Set_Name("sphere (mu=0.05)");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body("sphere");
        rigid_body->position=VECTOR_3D<T>(2,1,45);
        rigid_body->velocity=VECTOR_3D<T>(0,0,-10);
        rigid_body->angular_velocity=VECTOR_3D<T>(0,0,-10);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Set_Coefficient_Of_Friction(0.2);
        rigid_body->Set_Name("sphere (mu=0.05)");
        Append_Diffuse_Hints(color1);
    }
    else 
    {
        rigid_body=Initialize_Rigid_Body("sphere");
        rigid_body->position=VECTOR_3D<T>(-8,1,25);
        rigid_body->velocity=VECTOR_3D<T>(0,0,-10);
        rigid_body->Set_Coefficient_Of_Friction(0.05);
        rigid_body->Set_Name("sphere (mu=0.05)");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body("sphere");
        rigid_body->position=VECTOR_3D<T>(8,1,25);
        rigid_body->velocity=VECTOR_3D<T>(0,0,-10);
        rigid_body->Set_Coefficient_Of_Friction(0.02);
        rigid_body->Set_Name("sphere (mu=0.02)");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body("sphere");
        rigid_body->position=VECTOR_3D<T>(-2,10,0);
        rigid_body->angular_velocity=VECTOR_3D<T>(-10,0,0);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Set_Coefficient_Of_Friction(0.5);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Name("spinning sphere");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body("sphere");
        rigid_body->position=VECTOR_3D<T>(2,10,0);
        rigid_body->angular_velocity=VECTOR_3D<T>(-10,0,0);
        rigid_body->Update_Angular_Momentum();
        rigid_body->Set_Coefficient_Of_Friction(0.05);
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Name("low friction spinning sphere");
        Append_Diffuse_Hints(color1);
    }

    // PLANE -- no gravity
    rigid_body=Initialize_Rigid_Body("ground");
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    rigid_body->Set_Coefficient_Of_Restitution(1);
    Append_Diffuse_Hints(ground_color,false);
}

namespace PhysBAM
{
    template ROLLING_FRICTION_EXAMPLE<double>;
    template ROLLING_FRICTION_EXAMPLE<float>;
}
