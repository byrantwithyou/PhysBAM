//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Guendelman - January 4, 2003
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include "NONCONVEX_BOUNCE_EXAMPLE.h"

using namespace PhysBAM;

template<class T>
NONCONVEX_BOUNCE_EXAMPLE<T>::NONCONVEX_BOUNCE_EXAMPLE(int parameter) : RIGID_BODIES_EXAMPLE<T>(parameter)
{
    last_frame=240;
    output_directory="Nonconvex_Bounce_Example/output";
}

template<class T>
NONCONVEX_BOUNCE_EXAMPLE<T>::~NONCONVEX_BOUNCE_EXAMPLE()
{}

//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void NONCONVEX_BOUNCE_EXAMPLE<T>::
Initialize_Rigid_Bodies()
{
    VECTOR_3D<double> ground_color=RGB_COLORS<T>().Ground_Tan(),color1=RGB_COLORS<T>().Light_Steel_Blue();

    RIGID_BODY<TV> *rigid_body = 0;

    T scale = 20;
    T mass = 1;
    char bonefile[256];
    T height;
    T separation;

    if (parameter == 1) {
        strcpy(bonefile, "New_Bones/Left_Femur_2");
        height = 10;
        separation = 15;
    }
    else {
        strcpy(bonefile, "New_Bones/Pelvis_1"); 
        height = 5;
        separation = 3;
    }

    rigid_body=Initialize_Rigid_Body(bonefile, scale);
    rigid_body->Set_Mass(mass);
    rigid_body->position=VECTOR_3D<T>(-separation,height,0);
    if (parameter == 0) rigid_body->orientation=QUATERNION<T>(0.05,VECTOR_3D<T>(0,0,1))*QUATERNION<T>(-0.5,VECTOR_3D<T>(1,0,0))*QUATERNION<T>(pi,VECTOR_3D<T>(0,1,0));
    else rigid_body->orientation = QUATERNION<T>(pi/2,VECTOR_3D<T>(0,0,1));
    rigid_body->Set_Coefficient_Of_Restitution(1.0);
    rigid_body->Set_Name("bone (coeff 1)");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body(bonefile, scale);
    rigid_body->Set_Mass(mass);
    rigid_body->position=VECTOR_3D<T>(0,height,0);
    if (parameter == 0) rigid_body->orientation=QUATERNION<T>(0.02,VECTOR_3D<T>(0,0,1))*QUATERNION<T>(-0.5,VECTOR_3D<T>(1,0,0))*QUATERNION<T>(pi,VECTOR_3D<T>(0,1,0));
    else rigid_body->orientation = QUATERNION<T>(pi/2,VECTOR_3D<T>(0,0,1));
    rigid_body->Set_Coefficient_Of_Restitution(0.6);
    rigid_body->Set_Name("bone (coeff 0.6)");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body(bonefile, scale);
    rigid_body->Set_Mass(mass);
    rigid_body->position=VECTOR_3D<T>(separation,height,0);
    if (parameter == 0) rigid_body->orientation=QUATERNION<T>(0,VECTOR_3D<T>(0,0,1))*QUATERNION<T>(-0.5,VECTOR_3D<T>(1,0,0))*QUATERNION<T>(pi,VECTOR_3D<T>(0,1,0));
    else rigid_body->orientation = QUATERNION<T>(pi/2,VECTOR_3D<T>(0,0,1));
    rigid_body->Set_Coefficient_Of_Restitution(0.0);
    rigid_body->Set_Name("bone (coeff 0)");
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
    template NONCONVEX_BOUNCE_EXAMPLE<double>;
    template NONCONVEX_BOUNCE_EXAMPLE<float>;
}
