//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/KD_TREE.h>
#include <Core/Matrices/FRAME.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Basic_Geometry/CONE.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Basic_Geometry/TORUS.h>
#include <Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_DILATE.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INTERSECTION.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INVERT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_UNION.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_UTILITIES.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_SPHERE.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_VISCOSITY.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <fstream>
#include "STANDARD_TESTS_3D.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,3> >::
STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_BASE<TV>(stream_type_input,parse_args),
    foo_int1(0),foo_T1(0),foo_T2(0),foo_T3(0),foo_T4(0),foo_T5(0),
    use_foo_T1(false),use_foo_T2(false),use_foo_T3(false),use_foo_T4(false),use_foo_T5(false)
{
    parse_args.Add("-fooint1",&foo_int1,"int1","a interger");
    parse_args.Add("-fooT1",&foo_T1,&use_foo_T1,"T1","a scalar");
    parse_args.Add("-fooT2",&foo_T2,&use_foo_T2,"T2","a scalar");
    parse_args.Add("-fooT3",&foo_T3,&use_foo_T3,"T3","a scalar");
    parse_args.Add("-fooT4",&foo_T4,&use_foo_T4,"T4","a scalar");
    parse_args.Add("-fooT5",&foo_T5,&use_foo_T5,"T5","a scalar");
    parse_args.Parse();
    if(!this->override_output_directory) output_directory=LOG::sprintf("Test_3d_%i",test_number);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,3> >::
~STANDARD_TESTS()
{
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Write_Output_Files(const int frame)
{
    if(write_output_files) write_output_files(frame);
    BASE::Write_Output_Files(frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Read_Output_Files(const int frame)
{
    if(read_output_files) read_output_files(frame);
    BASE::Read_Output_Files(frame);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Initialize()
{
    switch(test_number)
    {
        case 1:{ // rotating sphere
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            SPHERE<TV> sphere(TV(.5,.5,.5)*m,.3*m);
            VECTOR<T,3> angular_velocity(TV(0.4,0,0)/s);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return angular_velocity.Cross(X-sphere.center);},
                [=](const TV&){return MATRIX<T,3>::Cross_Product_Matrix(angular_velocity);}
                ,density,particles_per_cell);
            T total_mass=particles.mass.Sum();
            TV total_momentum=particles.V.Weighted_Sum(particles.mass);
            TV dV=total_momentum/total_mass;
            particles.V-=dV;
            Add_Quasi_Pressure(1e3*unit_p*scale_E,7);
        } break;
        case 2:{ // Oscillating sphere
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            SPHERE<TV> sphere(TV(.5,.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,[=](const TV&){return MATRIX<T,3>();},
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,3>()+1.5);
            Add_Quasi_Pressure(1e3*unit_p*scale_E,7);
        } break;
        case 3:{ // Freefall sphere
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            SPHERE<TV> sphere(TV(.5,.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,[=](const TV&){return MATRIX<T,3>();},
                density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;
        case 7:{ // skew impact of two elastic spheres
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(),TV(30,30,30))*m,true);
            T density=5*unit_rho*scale_mass;
            SPHERE<TV> sphere1(TV(10,13,15)*m,2*m);
            Seed_Particles(sphere1,[=](const TV& X){return TV(0.75,0,0)*(m/s);},[=](const TV&){return MATRIX<T,3>();},density,particles_per_cell);
            SPHERE<TV> sphere2(TV(20,15,15)*m,2*m);
            Seed_Particles(sphere2,[=](const TV& X){return TV(-0.75,0,0)*(m/s);},[=](const TV&){return MATRIX<T,3>();},density,particles_per_cell);
            Add_Quasi_Pressure(31.685*unit_p*scale_E,7);
        } break;
        case 11:{ // skew impact of two elastic spheres with initial angular velocity
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(),TV(30,30,30))*m,true);
            T density=5*unit_rho*scale_mass;
            SPHERE<TV> sphere1(TV(10,13,15)*m,2*m);
            VECTOR<T,3> angular_velocity1(TV(0,0,foo_T1));
            Seed_Particles(sphere1,[=](const TV& X){return angular_velocity1.Cross(X-sphere1.center)+TV(0.75,0,0)*(m/s);},
                [=](const TV&){return MATRIX<T,3>::Cross_Product_Matrix(angular_velocity1);},density,particles_per_cell);
            SPHERE<TV> sphere2(TV(20,15,15)*m,2*m);
            VECTOR<T,3> angular_velocity2(TV(0,0,foo_T2));
            Seed_Particles(sphere2,[=](const TV& X){return angular_velocity2.Cross(X-sphere2.center)+TV(-0.75,0,0)*(m/s);},
                [=](const TV&){return MATRIX<T,3>::Cross_Product_Matrix(angular_velocity2);},density,particles_per_cell);
            Add_Quasi_Pressure(31.685*unit_p*scale_E,7);
        } break;
        case 30:{ // (fluid test) pool of water 
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            RANGE<TV> box(TV(0,0,0)*m,TV(1,0.25,1)*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,[=](const TV&){return MATRIX<T,3>();},density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;
        case 31:{ // (fluid test) circle drop 
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            SPHERE<TV> sphere(TV(.5,.75,.5)*m,.2*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,[=](const TV&){return MATRIX<T,3>();},density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-9.8,0));
        } break;
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
}
template class STANDARD_TESTS<VECTOR<float,3> >;
template class STANDARD_TESTS<VECTOR<double,3> >;
}
