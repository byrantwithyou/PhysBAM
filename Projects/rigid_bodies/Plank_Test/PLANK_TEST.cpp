//#####################################################################
// Copyright 2003, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Guendelman - January 2, 2003
//#####################################################################
#include "PLANK_TEST.h"

using namespace PhysBAM;

template<class T>
PLANK_TEST<T>::PLANK_TEST(int parameter) : RIGID_BODIES_EXAMPLE<T>(parameter)
{
    last_frame=480;
    output_directory="Plank_Test/output";

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
PLANK_TEST<T>::~PLANK_TEST()
{}

//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void PLANK_TEST<T>::
Initialize_Rigid_Bodies()
{
    VECTOR_3D<double> ground_color=RGB_COLORS<T>().Ground_Tan(),color1=RGB_COLORS<T>().Light_Steel_Blue();

    RIGID_BODY<TV> *rigid_body = 0;

    T baseboxsize=2;
    T dropboxsize=1.2;
    T plankscale=1;
    T dropheight=70;
    T smallboxsize=0.5;
    T smallboxmass=1;
    T offset = 0.05;
    const char *boxfile = (parameter == 0) ? "subdivided_box" : "box";
    const char *plankfile = (parameter == 0) ? "plank" : "unsubdivided_plank";

    T stack_epsilon = 0.3;
    T stack_mu = 0.5;

    rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
    rigid_body->position=VECTOR_3D<T>(-0.65, baseboxsize+4, 0);
    rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body->Set_Coefficient_Of_Friction(stack_mu);
    rigid_body->Set_Mass(smallboxmass);
    rigid_body->Set_Name("stack box 1a");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
    rigid_body->position=VECTOR_3D<T>(0.65, baseboxsize+5, 0);
    rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body->Set_Coefficient_Of_Friction(stack_mu);
    rigid_body->Set_Mass(smallboxmass);
    rigid_body->Set_Name("stack box 1b");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
    rigid_body->position=VECTOR_3D<T>(offset, baseboxsize+7, 0);
    rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body->Set_Coefficient_Of_Friction(stack_mu);
    rigid_body->Set_Mass(smallboxmass);
    rigid_body->Set_Name("stack box 2");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
    rigid_body->position=VECTOR_3D<T>(0, baseboxsize+9, offset);
    rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body->Set_Coefficient_Of_Friction(stack_mu);
    rigid_body->Set_Mass(smallboxmass);
    rigid_body->Set_Name("stack box 3");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
    rigid_body->position=VECTOR_3D<T>(-offset, baseboxsize+14, -offset);
    rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body->Set_Coefficient_Of_Friction(stack_mu);
    rigid_body->Set_Mass(smallboxmass);
    rigid_body->Set_Name("stack box 4");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
    rigid_body->position=VECTOR_3D<T>(0, baseboxsize+17, 0);
    rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body->Set_Coefficient_Of_Friction(stack_mu);
    rigid_body->Set_Mass(smallboxmass);
    rigid_body->Set_Name("stack box 5");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
    rigid_body->position=VECTOR_3D<T>(0, baseboxsize+19, 0);
    rigid_body->orientation=QUATERNION<T>(pi/4,VECTOR_3D<T>(0,1,0));
    rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body->Set_Coefficient_Of_Friction(stack_mu);
    rigid_body->Set_Mass(smallboxmass);
    rigid_body->Set_Name("stack box 6");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body(boxfile);
    rigid_body->position=VECTOR_3D<T>(baseboxsize+plankscale*5-dropboxsize,
                                   2*baseboxsize+0.5*plankscale+dropboxsize+dropheight,0);
    rigid_body->orientation=QUATERNION<T>(pi/2,VECTOR_3D<T>(0,1,0));
    rigid_body->Set_Coefficient_Of_Restitution(0.05);
    rigid_body->Set_Mass(100);
    rigid_body->Set_Name("drop box");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body(plankfile,plankscale);
    rigid_body->position=VECTOR_3D<T>(baseboxsize,2*baseboxsize+0.25,0);
    rigid_body->orientation=QUATERNION<T>(pi/2,VECTOR_3D<T>(0,1,0));
    rigid_body->Set_Coefficient_Of_Restitution(0.1);
    rigid_body->Set_Mass(10);
    rigid_body->Set_Name("plank");
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body(boxfile,baseboxsize);
    rigid_body->position=VECTOR_3D<T>(0,baseboxsize,0);
    rigid_body->Set_Coefficient_Of_Restitution(0.1);
    rigid_body->Set_Coefficient_Of_Friction(0.3);
    rigid_body->Set_Name("base box");
    rigid_body->is_static = true;
    Append_Diffuse_Hints(color1);

    rigid_body=Initialize_Rigid_Body("ground");
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    rigid_body->Set_Coefficient_Of_Restitution(1);
    Append_Diffuse_Hints(ground_color,false);
}

namespace PhysBAM
{
    template PLANK_TEST<double>;
    template PLANK_TEST<float>;
}
