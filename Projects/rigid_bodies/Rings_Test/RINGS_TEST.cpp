//#####################################################################
// Copyright 2002, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Guendelman - December 3, 2002
//#####################################################################
#include "RINGS_TEST.h"

using namespace PhysBAM;

template<class T>
RINGS_TEST<T>::RINGS_TEST(int parameter) : RIGID_BODIES_EXAMPLE<T>(parameter)
{
    last_frame=720;
    output_directory="Rings_Test/output";
    artificial_maximum_speed=40;
}

template<class T>
RINGS_TEST<T>::~RINGS_TEST()
{}

//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void RINGS_TEST<T>::
Initialize_Rigid_Bodies()
{
    VECTOR_3D<double> ground_color=RGB_COLORS<T>().Ground_Tan();
    VECTOR_3D<double> ring_color=RGB_COLORS<T>().Light_Steel_Blue();
    VECTOR_3D<double> pole_color=RGB_COLORS<T>().Gold();

    RIGID_BODY<TV> *rigid_body = 0;
    T mu = 0.6;
    T epsilon = 0.3;
    int i;

    if (parameter == 0)
    {
        std::cout << "Parameter 0: running 500 rings" << std::endl;
        CYLINDRICAL_RANDOM_PLACEMENT<T> random_placement(20, 300, VECTOR_3D<T>(0,30,0));
        random_placement.Set_Max_Orientation_Angle(0.2);
        random_placement.Set_Speed_Range(0, 1);
        random_placement.Set_Angular_Speed_Range(0, 1);
        Random_Scene_Generator("Rings_Test/ring_revolve", 500, ring_color, 12345, random_placement);
    }
    else
    {
        std::cout << "Parameter 1: running 1000 rings" << std::endl;
        CYLINDRICAL_RANDOM_PLACEMENT<T> random_placement(20, 600, VECTOR_3D<T>(0,30,0));
        random_placement.Set_Max_Orientation_Angle(0.2);
        random_placement.Set_Speed_Range(0, 1);
        random_placement.Set_Angular_Speed_Range(0, 1);
        Random_Scene_Generator("Rings_Test/ring_revolve", 1000, ring_color, 11111, random_placement);
    }

    for (i = 1; i <= rigid_body_parameters.list.rigid_bodies.m; i++)
    {
        rigid_body_parameters.list.rigid_bodies(i)->Set_Coefficient_Of_Restitution(epsilon);
        rigid_body_parameters.list.rigid_bodies(i)->Set_Coefficient_Of_Friction(mu);
        rigid_body_parameters.list.rigid_bodies(i)->Set_Mass(10);
        rigid_body_parameters.list.rigid_bodies(i)->triangulated_surface->Set_Desired_Particle_Partition_Size(3,1,3);
    }

    int poles = 5;

    for (i = 1; i <= poles; i++)
        for (int j = 1; j <= poles; j++)
        {
            // Poles
            rigid_body=Initialize_Rigid_Body("Rings_Test/medium_cylinder");
            char name[256];
            sprintf(name, "pole %d %d", i, j);
            rigid_body->Set_Name(name);
            rigid_body->is_static=true;
            rigid_body->position = VECTOR_3D<T>((i-(poles+1)/2.0)*7,10,(j-(poles+1)/2.0)*7);
            rigid_body->triangulated_surface->Set_Desired_Particle_Partition_Size(1,10,1);
            rigid_body->Set_Coefficient_Of_Friction(mu);
            Append_Diffuse_Hints(pole_color,false);
        }

    // PLANE -- no gravity
    rigid_body=Initialize_Rigid_Body("ground",2);
    rigid_body->Set_Coefficient_Of_Friction(mu);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    rigid_body->triangulated_surface->Set_Desired_Particle_Partition_Size(1,1,1);
    Append_Diffuse_Hints(ground_color,false);
}

namespace PhysBAM
{
    template RINGS_TEST<double>;
    template RINGS_TEST<float>;
}
