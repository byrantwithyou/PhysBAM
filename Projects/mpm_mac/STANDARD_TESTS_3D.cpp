//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INTERSECTION.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INVERT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_UNION.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <fstream>
#include "STANDARD_TESTS_3D.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,3> >::
STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_BASE<TV>(stream_type_input,parse_args)
{
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
// Function Initialize
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Initialize()
{
    switch(test_number)
    {
        case 1:{ // rotating sphere
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5,.5)*m,.3*m);
            VECTOR<T,3> angular_velocity(TV(0.4,0,0)/s);
            density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return angular_velocity.Cross(X-sphere.center);},
                [=](const TV&){return MATRIX<T,3>::Cross_Product_Matrix(angular_velocity);}
                ,density,particles_per_cell);
            T total_mass=particles.mass.Sum();
            TV total_momentum=particles.V.Weighted_Sum(particles.mass);
            TV dV=total_momentum/total_mass;
            particles.V-=dV;
        } break;
        case 2:{ // Oscillating sphere
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5,.5)*m,.3*m);
            density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,[=](const TV&){return MATRIX<T,3>();},
                density,particles_per_cell);
        } break;
        case 3:{ // Freefall sphere
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5,.5)*m,.3*m);
            density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,[=](const TV&){return MATRIX<T,3>();},
                density,particles_per_cell);
        } break;
        case 4:{ // Dam break
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            density=unit_rho*scale_mass;
            Seed_Particles(RANGE<TV>(TV(.0,.0,.0),TV(m/3,2*m/3,m)),0,0,density,particles_per_cell);
            gravity=TV(0,-1,0)*m/sqr(s);
        } break;
        case 100: this->Commandline_Analytic_test(); break;
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
    if(forced_collision_type!=-1)
        for(int i=0;i<collision_objects.m;i++)
            collision_objects(i)->type=(COLLISION_TYPE)forced_collision_type;
}
template class STANDARD_TESTS<VECTOR<float,3> >;
template class STANDARD_TESTS<VECTOR<double,3> >;
}
