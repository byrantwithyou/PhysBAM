//#####################################################################
// Copyright 2002, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Guendelman - December 3, 2002
//#####################################################################
#include <iostream>
#include "SLIDING_TEST.h"

using namespace PhysBAM;

template<class T>
SLIDING_TEST<T>::SLIDING_TEST(int parameter) : RIGID_BODIES_EXAMPLE<T>(parameter)
{
    last_frame=480;
    output_directory="Sliding_Test/output";

    cfl = 0.05; // for higher accuracy

    if (parameter == 0) {

    }
    else if (parameter == 1) {
//    no longer supported
//        integrate_velocity_before_collisions_one_body = 2;
    }
}

template<class T>
SLIDING_TEST<T>::~SLIDING_TEST()
{}

//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void SLIDING_TEST<T>::
Initialize_Rigid_Bodies()
{
    VECTOR_3D<double> ground_color=RGB_COLORS<T>().Ground_Tan(),color1=RGB_COLORS<T>().Light_Steel_Blue();

    RIGID_BODY<TV> *rigid_body = 0;

    char name[256];
    const char *boxfile = "box";
    T ground_angle_deg = 20;
    T box_size = 0.05;

    T epsilon1 = 1;
    T mu1 = 0.5;
    T initial_speed_down_ramp1 = 1.25;

    T epsilon2 = 1;
    T mu2 = 0.25;
    T initial_speed_down_ramp2 = 0;

    T ground_angle_rad = ground_angle_deg * pi / 180.0;
    std::cout << "Using angle " << ground_angle_deg << " deg (" << ground_angle_rad << ") (critical coeff=" << tan(ground_angle_rad) << std::endl;
    std::cout << "mu1 = " << mu1 << ", epsilon1 = " << epsilon1 << std::endl;
    std::cout << "mu2 = " << mu2 << ", epsilon2 = " << epsilon2 << std::endl;

    rigid_body=Initialize_Rigid_Body(boxfile,box_size);
    rigid_body->position=VECTOR_3D<T>(-box_size*sin(ground_angle_rad), box_size*cos(ground_angle_rad), 0);
    rigid_body->orientation = QUATERNION<T>(ground_angle_rad, VECTOR_3D<T>(0,0,1));
    rigid_body->velocity = rigid_body->orientation.Rotate(VECTOR_3D<T>(-initial_speed_down_ramp1,0,0));
    rigid_body->Set_Coefficient_Of_Friction(mu1);
    rigid_body->Set_Coefficient_Of_Restitution(epsilon1);
    sprintf(name, "box 1 epsilon=%.2f mu=%.2f", epsilon1, mu1); rigid_body->Set_Name(name);
    Append_Diffuse_Hints(color1);

    if (parameter == 0) {
        rigid_body=Initialize_Rigid_Body(boxfile,box_size);
        rigid_body->position=VECTOR_3D<T>(-box_size*sin(ground_angle_rad), box_size*cos(ground_angle_rad), 8*box_size);
        rigid_body->orientation = QUATERNION<T>(ground_angle_rad, VECTOR_3D<T>(0,0,1));
        rigid_body->velocity = rigid_body->orientation.Rotate(VECTOR_3D<T>(-initial_speed_down_ramp2,0,0));
        rigid_body->Set_Coefficient_Of_Friction(mu2);
        rigid_body->Set_Coefficient_Of_Restitution(epsilon2);
        sprintf(name, "box 2 epsilon=%.2f mu=%.2f", epsilon2, mu2); rigid_body->Set_Name(name);
        Append_Diffuse_Hints(color1);
    }
    else if (parameter == 1) {
        rigid_body=Initialize_Rigid_Body(boxfile,box_size);
        rigid_body->position=VECTOR_3D<T>(-box_size*sin(ground_angle_rad), box_size*cos(ground_angle_rad), 8*box_size);
        rigid_body->orientation = QUATERNION<T>(ground_angle_rad, VECTOR_3D<T>(0,0,1));
        rigid_body->velocity = rigid_body->orientation.Rotate(VECTOR_3D<T>(-initial_speed_down_ramp1,0,0));
        rigid_body->Set_Coefficient_Of_Friction(mu1);
        rigid_body->Set_Coefficient_Of_Restitution(epsilon1);
        sprintf(name, "box 2 velocity integrated before collisions"); rigid_body->Set_Name(name);
        Append_Diffuse_Hints(color1);
    }

    // PLANE -- no gravity
    rigid_body=Initialize_Rigid_Body("ground");
    rigid_body->Set_Name("ground");
    rigid_body->orientation = QUATERNION<T>(ground_angle_rad, VECTOR_3D<T>(0,0,1));
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    rigid_body->Set_Coefficient_Of_Friction(max(mu1,mu2));
    rigid_body->Set_Coefficient_Of_Restitution(max(epsilon1,epsilon2));
    Append_Diffuse_Hints(ground_color,false);
}

namespace PhysBAM
{
    template SLIDING_TEST<double>;
    template SLIDING_TEST<float>;
}
