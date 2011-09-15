//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Guendelman - May 2, 2003
//#####################################################################
#include "SEESAW_EXAMPLE.h"

using namespace PhysBAM;

template<class T>
SEESAW_EXAMPLE<T>::SEESAW_EXAMPLE(int parameter) : RIGID_BODIES_EXAMPLE<T>(parameter)
{
    last_frame=480;
    output_directory="Seesaw_Example/output";
    artificial_maximum_speed=20;

    if (parameter == 0) {

    }
    else if (parameter == 1) {
        use_triangle_hierarchy = true;
        use_edge_intersection = true;
        use_triangle_hierarchy_center_phi_test = true;
    }
}

template<class T>
SEESAW_EXAMPLE<T>::~SEESAW_EXAMPLE()
{}

//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void SEESAW_EXAMPLE<T>::
Initialize_Rigid_Bodies()
{
    VECTOR_3D<double> ground_color=RGB_COLORS<T>().Ground_Tan(),color1=RGB_COLORS<T>().Light_Steel_Blue();
    RIGID_BODY<TV> *rigid_body = 0;
    T mu = 0.2;
    bool fromrest = true;

    rigid_body=Initialize_Rigid_Body("subdivided_box",0.9);
    rigid_body->Set_Mass(20);
    if (fromrest)
        rigid_body->position=VECTOR_3D<T>(3.5,1.5+0.9,0);
    else
        rigid_body->position=VECTOR_3D<T>(3.5,20,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.05);
    rigid_body->Set_Coefficient_Of_Friction(mu);
    rigid_body->Set_Name("box1");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body("subdivided_box",0.5);
    rigid_body->Set_Mass(5);
    if (fromrest)
        rigid_body->position=VECTOR_3D<T>(-4,1.5+0.5,0);
    else
        rigid_body->position=VECTOR_3D<T>(-4,3,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.05);
    rigid_body->Set_Coefficient_Of_Friction(mu);
    rigid_body->Set_Name("box2");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body("plank");
    if (fromrest)
        rigid_body->position=VECTOR_3D<T>(0,1.25,0);
    else
        rigid_body->position=VECTOR_3D<T>(0,1.5,0);
    rigid_body->orientation=QUATERNION<T>(pi/2,VECTOR_3D<T>(0,1,0));
    rigid_body->Set_Coefficient_Of_Restitution(0.05);
    rigid_body->Set_Mass(10);
    rigid_body->Set_Name("plank");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body("Rings_Test/cylinder_revolve");
    rigid_body->position=VECTOR_3D<T>(0,0.5,0);
    rigid_body->orientation=QUATERNION<T>(pi/2,VECTOR_3D<T>(1,0,0));
    rigid_body->Set_Name("cylinder");
    rigid_body->is_static=true;
    rigid_body->Set_Coefficient_Of_Restitution(1);
    Append_Diffuse_Hints(color1);

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
    template SEESAW_EXAMPLE<double>;
    template SEESAW_EXAMPLE<float>;
}
