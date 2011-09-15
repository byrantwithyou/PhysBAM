//#####################################################################
// Copyright 2002,2003, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Guendelman - January 1, 2003
//#####################################################################
#include "STACKING_EXAMPLE.h"

using namespace PhysBAM;

template<class T>
STACKING_EXAMPLE<T>::STACKING_EXAMPLE(int parameter) : RIGID_BODIES_EXAMPLE<T>(parameter)
{
    last_frame=1440;
    output_directory="Stacking_Example/output";

    artificial_maximum_speed=5;

    max_rotation_per_time_step = 0.4*pi;
    max_linear_movement_fraction_per_time_step = 0.4;
}

template<class T>
STACKING_EXAMPLE<T>::~STACKING_EXAMPLE()
{}

//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void STACKING_EXAMPLE<T>::
Initialize_Rigid_Bodies()
{
    RIGID_BODY<TV> *rigid_body=0;

    VECTOR_3D<double> ground_color=RGB_COLORS<T>().Ground_Tan();
    VECTOR_3D<double> color1=RGB_COLORS<T>().Light_Steel_Blue();
    VECTOR_3D<double> color2=RGB_COLORS<T>().Red();

    ARRAY<std::string>filenames;
    int i;
    for (i = 1; i <= 300; i++) filenames.Append("New_Bones/Cranium_1");
    for (i = 1; i <= 150; i++) filenames.Append("New_Bones/Pelvis_1");
    for (i = 1; i <= 50; i++) filenames.Append("New_Bones/Left_Femur_2");

    T mu = 0.4;
    T epsilon = 0.3;
    T rolling_friction = 0.1;

    CYLINDRICAL_RANDOM_PLACEMENT<T> random_placement(1, 40, VECTOR_3D<T>(0,3,0));
    random_placement.Set_Fixed_Scale(3);
    Random_Scene_Generator(filenames, color1, 12345, random_placement);
    T mass_scale = 10;
    for (i = 1; i <= rigid_body_parameters.list.rigid_bodies.m; i++)
    {
        rigid_body_parameters.list.rigid_bodies(i)->Set_Mass(rigid_body_parameters.list.rigid_bodies(i)->mass * mass_scale);
        rigid_body_parameters.list.rigid_bodies(i)->Set_Coefficient_Of_Restitution(epsilon);
        rigid_body_parameters.list.rigid_bodies(i)->Set_Coefficient_Of_Friction(mu);
        rigid_body_parameters.list.rigid_bodies(i)->Set_Coefficient_Of_Rolling_Friction(rolling_friction);

        if (rigid_body_parameters.list.rigid_bodies(i)->name.find("Cranium"))
            rigid_body_parameters.list.rigid_bodies(i)->triangulated_surface->Set_Desired_Particle_Partition_Size(2,2,2);
        else if (rigid_body_parameters.list.rigid_bodies(i)->name.find("Left_Femur"))
            rigid_body_parameters.list.rigid_bodies(i)->triangulated_surface->Set_Desired_Particle_Partition_Size(1,4,1); 
        else if (rigid_body_parameters.list.rigid_bodies(i)->name.find("Pelvis"))
            rigid_body_parameters.list.rigid_bodies(i)->triangulated_surface->Set_Desired_Particle_Partition_Size(3,3,2);
    }

    // Funnel
//    rigid_body=Initialize_Rigid_Body("funnel_thicker_revolve_new",0.4);
//    rigid_body->position=VECTOR_3D<double>(0,2,0);
    rigid_body=Initialize_Rigid_Body("funnel_thicker_revolve",0.4);
    rigid_body->position=VECTOR_3D<T>(0,2,0);
    rigid_body->Set_Name("funnel");
    rigid_body->is_static=true;
    rigid_body->triangulated_surface->Set_Desired_Particle_Partition_Size(3,3,3);
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);
    Append_Diffuse_Hints(color2);

    T yy = 0.2;
    T ss = 0.25;
    T xx = ss*(5-0.25);

    rigid_body=Initialize_Rigid_Body("plank",ss);
    rigid_body->position=VECTOR_3D<T>(xx,yy,0);
    rigid_body->orientation=QUATERNION<T>(pi/2,VECTOR_3D<T>(0,0,1));
    rigid_body->Set_Name("plank 1");
    rigid_body->is_static=true;
    rigid_body->triangulated_surface->Set_Desired_Particle_Partition_Size(1,10,1);
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);
    Append_Diffuse_Hints(color2);

    rigid_body=Initialize_Rigid_Body("plank",ss);
    rigid_body->position=VECTOR_3D<T>(0,yy,xx);
    rigid_body->orientation=QUATERNION<T>(pi/2,VECTOR_3D<T>(0,0,1));
    rigid_body->orientation=QUATERNION<T>(pi/2,VECTOR_3D<T>(0,1,0))*rigid_body->orientation;
    rigid_body->Set_Name("plank 2");
    rigid_body->is_static=true;
    rigid_body->triangulated_surface->Set_Desired_Particle_Partition_Size(1,10,1);
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);
    Append_Diffuse_Hints(color2);

    rigid_body=Initialize_Rigid_Body("plank",ss);
    rigid_body->position=VECTOR_3D<T>(-xx,yy,0);
    rigid_body->orientation=QUATERNION<T>(pi/2,VECTOR_3D<T>(0,0,1));
    rigid_body->orientation=QUATERNION<T>(pi,VECTOR_3D<T>(0,1,0))*rigid_body->orientation;
    rigid_body->Set_Name("plank 3");
    rigid_body->is_static=true;
    rigid_body->triangulated_surface->Set_Desired_Particle_Partition_Size(1,10,1);
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);
    Append_Diffuse_Hints(color2);

    rigid_body=Initialize_Rigid_Body("plank",ss);
    rigid_body->position=VECTOR_3D<T>(0,yy,-xx);
    rigid_body->orientation=QUATERNION<T>(pi/2,VECTOR_3D<T>(0,0,1));
    rigid_body->orientation=QUATERNION<T>(3.0*pi/2,VECTOR_3D<T>(0,1,0))*rigid_body->orientation;
    rigid_body->Set_Name("plank 4");
    rigid_body->is_static=true;
    rigid_body->triangulated_surface->Set_Desired_Particle_Partition_Size(1,10,1);
    rigid_body->Set_Coefficient_Of_Rolling_Friction(0);
    Append_Diffuse_Hints(color2);

    // Ground
    rigid_body=Initialize_Rigid_Body("ground");
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    rigid_body->triangulated_surface->Set_Desired_Particle_Partition_Size(1,1,1);
    rigid_body->Set_Coefficient_Of_Friction(mu);
    rigid_body->Set_Coefficient_Of_Rolling_Friction(rolling_friction);
    Append_Diffuse_Hints(ground_color, false);
}

namespace PhysBAM
{
    template STACKING_EXAMPLE<double>;
    template STACKING_EXAMPLE<float>;
}
