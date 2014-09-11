//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_2D
//#####################################################################
//    1. What is it
//#####################################################################
#ifndef __STANDARD_TESTS_2D__
#define __STANDARD_TESTS_2D__

#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Implicit_Objects_Uniform/SMOOTH_LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/BEZIER_SPLINE.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Bindings/LINEAR_BINDING.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_PENALTY.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Forces/BEZIER_C2_FORCE.h>
#include <Deformables/Forces/BEZIER_CURVATURE_FORCE.h>
#include <Deformables/Forces/ELASTIC_ETHER_DRAG.h>
#include <Deformables/Forces/RALEIGH_DAMPING_FORCE.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_EVOLUTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_OBJECTIVE.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <fstream>
#include "STANDARD_TESTS_BASE.h"
namespace PhysBAM{

template<class TV> class STANDARD_TESTS;
template<class T_input>
class STANDARD_TESTS<VECTOR<T_input,2> >:public STANDARD_TESTS_BASE<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    typedef STANDARD_TESTS_BASE<TV> BASE;
    using BASE::solid_body_collection;using BASE::stream_type;using BASE::solids_evolution;using BASE::parse_args;
    using BASE::test_number;using BASE::data_directory;using BASE::m;using BASE::s;using BASE::kg;
    using BASE::tests;using BASE::density;using BASE::Get_Initial_Data_After;using BASE::use_penalty_self_collisions;
    using BASE::Initialize_Bodies_After;using BASE::resolution;using BASE::stiffness_multiplier;
    using BASE::curvature_stiffness_multiplier;using BASE::point_curves;using BASE::kinematic_points;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input)
        :BASE(stream_type_input)
    {
    }

    virtual ~STANDARD_TESTS()
    {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
//    parse_args->Add("-opt",&here,"meaning");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    bool automatically_add_to_collision_structures=true;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    switch(test_number){
        case 1:case 2:{
            if(resolution==0) resolution=10;
            ARRAY<TV> X(resolution+1);
            for(int i=0;i<=resolution;i++){
                T t=(T)i/resolution;
                X(i)=TV(t,sin(t*2*pi));}
            BEZIER_SPLINE<TV,3>* spline=BEZIER_SPLINE<TV,3>::Create(particles);
            Smooth_Fit<TV>(*spline,X);
            deformable_body_collection.Add_Structure(spline);
            deformable_body_collection.Add_Structure(Create_Segmented_Curve(*spline,true));
            deformable_body_collection.Add_Force(new BEZIER_CURVATURE_FORCE<TV>(particles,*spline,curvature_stiffness_multiplier,stiffness_multiplier));
            for(int i=0;i<spline->control_points.m;i++)
                particles.mass.Subset(spline->control_points(i))+=(T).25*density/spline->control_points.m;
            kinematic_points.Append(0);
            kinematic_points.Append(resolution*3);
            point_curves.Resize(2);
            point_curves(0).Add_Control_Point(0,particles.X(0));
            point_curves(0).Add_Control_Point(1,particles.X(0));
            point_curves(1).Add_Control_Point(0,particles.X(resolution*3));
            point_curves(1).Add_Control_Point(1,particles.X(resolution*3)+TV(test_number==1?1:-.9,0));
            for(int i=4;i<particles.X.m;i+=3)
                deformable_body_collection.binding_list.Add_Binding(
                    new LINEAR_BINDING<TV,2>(particles,i,VECTOR<int,2>(i-2,i-1),VECTOR<T,2>(-1,2)));
            particles.mass(0)=FLT_MAX;
            particles.mass(resolution*3)=FLT_MAX;
            break;}
        default:
            LOG::cerr<<"Initial Data: Unrecognized test number "<<test_number<<std::endl;exit(1);}

    Get_Initial_Data_After(automatically_add_to_collision_structures);

    switch(test_number){
        case 1:
        case 2:
            break;
        default:
            LOG::cerr<<"Missing bodies implementation for test number "<<test_number<<std::endl;exit(1);}

    Initialize_Bodies_After();
}
};
}
#endif
