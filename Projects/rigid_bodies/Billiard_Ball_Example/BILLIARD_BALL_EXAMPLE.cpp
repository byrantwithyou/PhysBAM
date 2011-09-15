//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Guendelman - January 4, 2003
//#####################################################################
#include "BILLIARD_BALL_EXAMPLE.h"

using namespace PhysBAM;

template<class T>
BILLIARD_BALL_EXAMPLE<T>::BILLIARD_BALL_EXAMPLE(int parameter) : RIGID_BODIES_EXAMPLE<T>(parameter)
{
    last_frame=240;
    output_directory="Billiard_Ball_Example/output";
}

template<class T>
BILLIARD_BALL_EXAMPLE<T>::~BILLIARD_BALL_EXAMPLE()
{}

//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void BILLIARD_BALL_EXAMPLE<T>::
Initialize_Rigid_Bodies()
{
    VECTOR_3D<double> ground_color=RGB_COLORS<T>().Ground_Tan(),color1=RGB_COLORS<T>().Light_Steel_Blue();

    RIGID_BODY<TV> *rigid_body = 0;
    T mu = 0.3;

    int number = 5;
    for (int k = 1; k <= number; k++) {
        rigid_body=Initialize_Rigid_Body("sphere");
        rigid_body->position=VECTOR_3D<T>(2*(k-(number+1)/2),1,0);
        char name[256]; sprintf(name, "sphere %d", k); 
        rigid_body->Set_Name(name);
        rigid_body->Set_Coefficient_Of_Restitution(0.9);
        rigid_body->Set_Coefficient_Of_Friction(mu);
        Append_Diffuse_Hints(color1);
    }

    rigid_body=Initialize_Rigid_Body("sphere");
    rigid_body->position=VECTOR_3D<T>(-30,1,0);
    rigid_body->velocity=VECTOR_3D<T>(20,0,0);
    rigid_body->Set_Name("sphere");
    rigid_body->Set_Coefficient_Of_Restitution(0.9);
    rigid_body->Set_Coefficient_Of_Friction(0.5);
    Append_Diffuse_Hints(color1);

    // PLANE -- no gravity
    rigid_body=Initialize_Rigid_Body("ground");
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->Set_Coefficient_Of_Friction(1);
    rigid_body->Set_Coefficient_Of_Rolling_Friction(1);
    Append_Diffuse_Hints(ground_color,false);
}

namespace PhysBAM
{
    template BILLIARD_BALL_EXAMPLE<double>;
    template BILLIARD_BALL_EXAMPLE<float>;
}
