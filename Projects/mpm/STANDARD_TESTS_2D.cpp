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
        case 2:{ // oscillating circle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV(0.1,0);},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>()+1.5);
            particles.B.Fill(MATRIX<T,2>(1,2,3,10));
            Add_Fixed_Corotated(1e3*scale_E,0.3);
        } break;
        case 3:{ // freefall circle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            Add_Gravity(TV(0,-9.8));
        } break;
        case 4:{ // colliding of two rings
            grid.Initialize(TV_INT(48,48),RANGE<TV>(TV(),TV(0.48,0.48)),true);
            ARRAY<SPHERE<TV> > spheres; ARRAY<TV> v0; ARRAY<T> r;
            spheres.Append(SPHERE<TV>(TV(0.1,0.24),0.04)); spheres.Append(SPHERE<TV>(TV(0.4,0.24),0.04));
            v0.Append(TV(50,0)); v0.Append(TV(-50,0));
            r.Append(0.03); r.Append(0.03);
            for(int s=0;s<spheres.m;s++){
                SPHERE<TV>& sphere=spheres(s);
                T density=1010*scale_mass;
                GRID<TV> sg(TV_INT(52,52),sphere.Bounding_Box().Thickened(r(s)*(T)1.5));
                int last=particles.number;
                Seed_Particles(sphere,[=](const TV& X){return v0(s);},[=](const TV&){return MATRIX<T,2>();},density,sg);
                for(int k=last;k<particles.number;k++){
                    if((particles.X(k)-sphere.center).Magnitude_Squared()<sqr(r(s))){
                        particles.deletion_list.Append(k);}}}
            particles.Delete_Elements_On_Deletion_List();
            Add_Neo_Hookean(0.073*1e9*scale_E,0.4);
        } break;
        case 5:{ // rebound of an elastic cylinder
            T dx=0.5; grid.Initialize(TV_INT(31+8,11+8),RANGE<TV>(TV(0-dx/2-dx*4,0-dx/2-dx*4),TV(15+dx/2+dx*4,5+dx/2+dx*4)),true);
            MPM_COLLISION_OBJECT<TV>* left=new MPM_COLLISION_OBJECT<TV>(new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(RANGE<TV>(TV(-5,-5),TV(0+dx/2,15))),
                new MPM_COLLISION_OBJECT_STATIC_PATH<TV>(FRAME<TV>()),false,0);
            MPM_COLLISION_OBJECT<TV>* right=new MPM_COLLISION_OBJECT<TV>(new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(RANGE<TV>(TV(15-dx/2,-5),TV(20,15))),
                new MPM_COLLISION_OBJECT_STATIC_PATH<TV>(FRAME<TV>()),false,0);
            collision_objects.Append(left);
            collision_objects.Append(right);
            SPHERE<TV> sphere(TV(2.5,2.5),1.5);
            T density=4*scale_mass;
            GRID<TV> sg(TV_INT(12,12),RANGE<TV>(TV(2.5-dx/2*5.5,2.5-dx/2*5.5),TV(2.5+dx/2*5.5,2.5+dx/2*5.5)));
            Seed_Particles(sphere,[=](const TV& X){return TV(0.5,0);},[=](const TV&){return MATRIX<T,2>();},
                density,sg);
            Add_Neo_Hookean(85.5*scale_E,0.425); //solve({E/(2*(1+r))=30,E*r/((1+r)*(1-2*r))=170},{E,r});
        } break;
        case 6:{ // skew impact of two elastic cylinders
            T dx=1; grid.Initialize(TV_INT(21+8,13+8),RANGE<TV>(TV(0-dx/2-dx*4,0-dx/2-dx*4),TV(20+dx/2+dx*4,12+dx/2+dx*4)),true);
            T density=5*scale_mass;
            SPHERE<TV> sphere1(TV(3,3),2);
            GRID<TV> sg1(TV_INT(12,12),RANGE<TV>(TV(2.5-dx/2*5.5,2.5-dx/2*5.5),TV(2.5+dx/2*5.5,2.5+dx/2*5.5)));
            Seed_Particles(sphere1,[=](const TV& X){return TV(0.75,0);},[=](const TV&){return MATRIX<T,2>();},density,sg1);
            SPHERE<TV> sphere2(TV(16,5),2);
            GRID<TV> sg2(TV_INT(12,12),RANGE<TV>(TV(2.5-dx/2*5+13.5,2.5-dx/2*5.5+2),TV(2.5+dx/2*5.5+13,2.5+dx/2*5.5+2)));
            Seed_Particles(sphere2,[=](const TV& X){return TV(-0.75,0);},[=](const TV&){return MATRIX<T,2>();},density,sg2);
            Add_Neo_Hookean(31.685*scale_E,0.44022); //solve({E/(2*(1+r))=11,E*r/((1+r)*(1-2*r))=81},{E,r});
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
