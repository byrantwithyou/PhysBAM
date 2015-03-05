//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/FRAME.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT_STATIC_PATH.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include "STANDARD_TESTS_2D.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,2> >::
STANDARD_TESTS(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_BASE<TV>(stream_type,parse_args)
{
    parse_args.Parse();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,2> >::
~STANDARD_TESTS()
{
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,2> >::
Write_Output_Files(const int frame)
{
    BASE::Write_Output_Files(frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,2> >::
Read_Output_Files(const int frame)
{
    BASE::Read_Output_Files(frame);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,2> >::
Initialize()
{
    switch(test_number)
    {
        case 1:{ // rotating circle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            VECTOR<T,1> angular_velocity(0.4);
            T density=2*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return angular_velocity.Cross(X-sphere.center);},
                [=](const TV&){return MATRIX<T,2>::Cross_Product_Matrix(angular_velocity);}
                ,density,particles_per_cell);
            T total_mass=particles.mass.Sum();
            TV total_momentum=particles.V.Weighted_Sum(particles.mass);
            TV dV=total_momentum/total_mass;
            particles.V-=dV;
            Add_Fixed_Corotated(1e3*scale_E,0.3);
        } break;
        case 2:{ // Oscillating circle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>()+1.5);
            Add_Fixed_Corotated(1e3*scale_E,0.3);
        } break;
        case 3:{ // Freefall circle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            Add_Gravity(TV(0,-9.8));
        } break;
        case 4:{ // Colliding spheres/rings
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            ARRAY<SPHERE<TV> > spheres;
            ARRAY<T> inner_radius;
            ARRAY<TV> V0;
            spheres.Append(SPHERE<TV>(TV(0.3,0.5),0.1));inner_radius.Append(0.05);V0.Append(TV(1,0));
            spheres.Append(SPHERE<TV>(TV(0.7,0.5),0.15));inner_radius.Append(0.03);V0.Append(TV(-0.5,0));
            for(int i=0;i<spheres.m;i++){
                T density=2*scale_mass;
                int last_m=particles.number;
                Seed_Particles(spheres(i),[=](const TV& X){return V0(i);},[=](const TV&){return MATRIX<T,2>();},
                    density,particles_per_cell);
                for(int k=last_m;k<particles.number;k++){
                    if((particles.X(k)-spheres(i).center).Magnitude_Squared()<sqr(inner_radius(i))){
                        particles.deletion_list.Append(k);}}}
            particles.Delete_Elements_On_Deletion_List();
            Add_Fixed_Corotated(1e3*scale_E,0.3);
        } break;
        case 5:{ // Dropping sphere to ground
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.2);
            T density=2*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            Add_Gravity(TV(0,-9.8));
            Add_Fixed_Corotated(1e3*scale_E,0.3);
            MPM_COLLISION_OBJECT<TV>* ground=new MPM_COLLISION_OBJECT<TV>(
                new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(RANGE<TV>(TV(),TV(1,.1))),
                new MPM_COLLISION_OBJECT_STATIC_PATH<TV>(FRAME<TV>()),true,10);
            collision_objects.Append(ground);
        } break;
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
}
//#####################################################################
// Function Begin_Frame
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,2> >::
Begin_Frame(const int frame)
{
}
//#####################################################################
// Function End_Frame
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,2> >::
End_Frame(const int frame)
{
}
//#####################################################################
// Function Begin_Time_Step
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,2> >::
Begin_Time_Step(const T time)
{
}
//#####################################################################
// Function End_Time_Step
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,2> >::
End_Time_Step(const T time)
{
}
template class STANDARD_TESTS<VECTOR<float,2> >;
template class STANDARD_TESTS<VECTOR<double,2> >;
}
