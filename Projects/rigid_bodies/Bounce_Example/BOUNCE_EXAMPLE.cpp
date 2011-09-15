//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Guendelman - January 2, 2003
//#####################################################################
#include "BOUNCE_EXAMPLE.h"

using namespace PhysBAM;

template<class T>
BOUNCE_EXAMPLE<T>::BOUNCE_EXAMPLE(int parameter) : RIGID_BODIES_EXAMPLE<T>(parameter)
{
    last_frame=240;
    output_directory="Bounce_Example/output";
}

template<class T>
BOUNCE_EXAMPLE<T>::~BOUNCE_EXAMPLE()
{}

//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void BOUNCE_EXAMPLE<T>::
Initialize_Rigid_Bodies()
{
    VECTOR_3D<double> ground_color=RGB_COLORS<T>().Ground_Tan(),color1=RGB_COLORS<T>().Light_Steel_Blue();

    RIGID_BODY<TV> *rigid_body = 0;

    // SPHERE
    rigid_body=Initialize_Rigid_Body("sphere");
    rigid_body->position=VECTOR_3D<T>(-3,5,0);
    rigid_body->Set_Coefficient_Of_Restitution(1.0);
    rigid_body->Set_Name("sphere (coeff 1)");
    Append_Diffuse_Hints(color1);

    // SPHERE
    rigid_body=Initialize_Rigid_Body("sphere");
    rigid_body->position=VECTOR_3D<T>(0,5,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.5);
    rigid_body->Set_Name("sphere (coeff 0.5)");
    Append_Diffuse_Hints(color1);

    // SPHERE
    rigid_body=Initialize_Rigid_Body("sphere");
    rigid_body->position=VECTOR_3D<T>(3,5,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.0);
    rigid_body->Set_Name("sphere (coeff 0)");
    Append_Diffuse_Hints(color1);

    // PLANE -- no gravity
    rigid_body=Initialize_Rigid_Body("ground");
    rigid_body->Set_Coefficient_Of_Restitution(1.0);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    Append_Diffuse_Hints(ground_color,false);
}

namespace PhysBAM
{
    template BOUNCE_EXAMPLE<double>;
    template BOUNCE_EXAMPLE<float>;
}
