//#####################################################################
// Copyright 2015, Lin Huang, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INTERSECTION.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INVERT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_UNION.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <fstream>
#include "STANDARD_TESTS_2D.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,2> >::
STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_BASE<TV>(stream_type_input,parse_args)
{
    parse_args.Parse();
    if(!this->override_output_directory) output_directory=LOG::sprintf("Test_%i",test_number);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,2> >::
~STANDARD_TESTS()
{
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,2> >::
Initialize()
{
    if(bc_periodic)
        bc_type.Fill(BC_PERIODIC);
    switch(test_number)
    {
        case 1:{ // stationary circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Set_Phases({density});
            Seed_Particles(sphere,0,0,density,particles_per_cell);
        } break;
        case 2:
        case 15:{ // translating circle
            if(test_number==2) bc_type.Fill(BC_INVALID);
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Set_Phases({density});
            Seed_Particles(sphere,[=](const TV& X){return TV(m/s,0);},0,density,particles_per_cell);
        } break;
        case 3:{ // rotating circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            VECTOR<T,1> angular_velocity(0.4/s);
            T density=2*unit_rho*scale_mass;
            Set_Phases({density});
            Seed_Particles(sphere,[=](const TV& X){return angular_velocity.Cross(X-sphere.center);},
                [=](const TV&){return MATRIX<T,2>::Cross_Product_Matrix(angular_velocity);},
                density,particles_per_cell);
            T total_mass=particles.mass.Sum();
            TV total_momentum=particles.V.Weighted_Sum(particles.mass);
            TV dV=total_momentum/total_mass;
            particles.V-=dV;
        } break;
        case 4:{ // freefall circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Set_Phases({density});
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
        } break;
        case 5:{ // stationary pool
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            Set_Phases({density});
            Seed_Particles(RANGE<TV>(TV(),TV(m,.5*m)),0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
        } break;
        case 6:{ // freefall rectangle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            Set_Phases({density});
            Seed_Particles(RANGE<TV>(TV(.2*m,.2*m),TV(.5*m,.8*m)),0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
        } break;
        case 7:{ // stationary circles in two phases
            T density=2*unit_rho*scale_mass;
            Set_Phases({density,density});
            particles.Store_Phase(true);
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere0(TV(.3,.3)*m,.1*m);
            Seed_Particles(sphere0,0,0,density,particles_per_cell);
            int n=particles.phase.m;
            SPHERE<TV> sphere1(TV(.6,.6)*m,.1*m);
            Seed_Particles(sphere1,0,0,density,particles_per_cell);
            particles.phase.Array_View(n,particles.phase.m-n).Fill(PHASE_ID(1));
        } break;
        case 8:{ // concave shape
            this->use_phi=true;
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            Set_Phases({density});
            SPHERE<TV> sphere0(TV(.3,.3)*m,.1*m);
            SPHERE<TV> sphere1(TV(.44,.44)*m,.1*m);
            auto shape=Unite(Make_IO(sphere0),Make_IO(sphere1));
            Seed_Particles(*shape,0,0,density,particles_per_cell);
            delete shape;
        } break;
        case 9: // freefall circles with different phases
        case 10:{ // freefall circles with same phase
            T density=2*unit_rho*scale_mass;
            if(test_number==9) Set_Phases({density,density});
            else Set_Phases({density});
            particles.Store_Phase(true);
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            gravity=TV(0,-1)*m/sqr(s);
            SPHERE<TV> sphere0(TV(.3,.5)*m,.1*m);
            SPHERE<TV> sphere1(TV(.7,.5)*m,.1*m);
            Seed_Particles(sphere0,0,0,density,particles_per_cell);
            int n=particles.phase.m;
            Seed_Particles(sphere1,0,0,density,particles_per_cell);
            if(test_number==9)
                particles.phase.Array_View(n,particles.phase.m-n).Fill(PHASE_ID(1));
            // a wall in the middle preventing the two circles from touching
            RANGE<TV> wall(TV(.45,0)*m,TV(.55,1)*m);
            Add_Collision_Object(wall,COLLISION_TYPE::slip,0);
        } break;
        case 11:
        case 12:{ // free fall circle with curved boundary
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            if(test_number==12) sphere.radius=.4*m;
            T density=2*unit_rho*scale_mass;
            Set_Phases({density});
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
            // add a circle collision object
            SPHERE<TV> sph(TV(.5,.5)*m,.4*m);
            Add_Collision_Object(Invert(Make_IO(sph)),COLLISION_TYPE::slip,0,0,0);
        } break;
        case 13:{ // half filled stationary pool with curved boundary
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            Set_Phases({density});
            SPHERE<TV> sphere(TV(.5,.5)*m,.4*m);
            RANGE<TV> box(TV(),TV(m,.5*m));
            auto seed_space=Intersect(Make_IO(sphere),Make_IO(box));
            Seed_Particles(*seed_space,0,0,density,particles_per_cell);
            delete seed_space;
            gravity=TV(0,-1)*m/sqr(s);
            // add a circle collision object
            SPHERE<TV> sph(TV(.5,.5)*m,.4*m);
            Add_Collision_Object(Invert(Make_IO(sph)),COLLISION_TYPE::slip,0,0,0);
        } break;
        case 14:{ // rotating circle (periodic test)
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            VECTOR<T,1> angular_velocity(0.4/s);
            T density=2*unit_rho*scale_mass;
            Set_Phases({density});
            for(int i=0;i<2;i++)
                for(int j=0;j<2;j++){
                    TV c(i*m,j*m);
                    auto io=Intersect(Make_IO(SPHERE<TV>(c,.3*m)),Make_IO(grid.domain));
                    Seed_Particles(*io,
                        [=](const TV& X){return angular_velocity.Cross(X-c);},
                        [=](const TV&){return MATRIX<T,2>::Cross_Product_Matrix(angular_velocity);},
                        density,particles_per_cell);
                    delete io;}
            T total_mass=particles.mass.Sum();
            TV total_momentum=particles.V.Weighted_Sum(particles.mass);
            TV dV=total_momentum/total_mass;
            particles.V-=dV;
            bc_type.Fill(BC_PERIODIC);
        } break;
        case 16:{ // Initialize velocity field as Taylor-Green Vortex
            // To specify parameters, use: -I a
            int a=extra_int.m>=1?extra_int(0):2;
            Set_Grid(RANGE<TV>::Unit_Box()*pi*m);
            T density=2*unit_rho*scale_mass;
            Set_Phases({density});
            auto V_func=[=](const TV& X){return TV(-sin(a*X.x)*cos(a*X.y),cos(a*X.x)*sin(a*X.y));};
            auto dV_func=[=](const TV& X)
                {
                    T c=cos(a*X.x)*cos(a*X.y),s=sin(a*X.x)*sin(a*X.y);
                    return MATRIX<T,2>(-c,-s,s,c)*a;
                };
            Seed_Particles(grid.domain,V_func,dV_func,density,particles_per_cell);
        } break;
        case 17:
        case 19:{ // stationary pool with two phases
            T water_density=1000*unit_rho*scale_mass;
            T air_density=unit_rho*scale_mass;
            Set_Phases({water_density,air_density});
            particles.Store_Phase(true);
            this->use_massless_particles=true;
            this->use_phi=true;
            this->use_multiphase_projection=true;
            if(test_number==19){
                this->ghost=4;
                this->use_bump=true;}
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            Seed_Particles(RANGE<TV>(TV(),TV(m,.5*m)),0,0,air_density,particles_per_cell);
            int n=particles.phase.m;
            Seed_Particles(RANGE<TV>(TV(0,.5*m),TV(m,m)),0,0,air_density,particles_per_cell);
            particles.phase.Array_View(n,particles.phase.m-n).Fill(PHASE_ID(1));
            gravity=TV(0,-1)*m/sqr(s);
        } break;
        case 18:{ // full of fluid with random inital velocities
            T a=extra_T.m>=1?extra_T(0):-1;
            T b=extra_T.m>=2?extra_T(1):1;
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            Set_Phases({density});
            RANDOM_NUMBERS<T> rand;
            rand.Set_Seed();
            auto V_func=[&](const TV&){
                return TV(rand.Get_Uniform_Number(a,b),rand.Get_Uniform_Number(a,b));};
            Seed_Particles(grid.domain,V_func,0,density,particles_per_cell);
        } break;
        case 20:{ // stationary sphere in air
            T water_density=1000*unit_rho*scale_mass;
            T air_density=unit_rho*scale_mass;
            Set_Phases({water_density,air_density});
            particles.Store_Phase(true);
            this->use_massless_particles=true;
            this->use_phi=true;
            this->use_multiphase_projection=true;
            this->ghost=4;
            this->use_bump=true;
            RANGE<TV> box=RANGE<TV>::Unit_Box()*m;
            Set_Grid(box);
            SPHERE<TV> sphere0(TV(.3,.3)*m,.1*m);
            SPHERE<TV> sphere1(TV(.43,.43)*m,.1*m);
            auto shape=Unite(Make_IO(sphere0),Make_IO(sphere1));
            Seed_Particles(*shape,0,0,water_density,particles_per_cell);
            int n=particles.phase.m;
            box.Scale_About_Center(0.9);
            auto cshape=Intersect(Make_IO(box),Invert(shape));
            Seed_Particles(*cshape,0,0,air_density,particles_per_cell);
            particles.phase.Array_View(n,particles.phase.m-n).Fill(PHASE_ID(1));
            delete cshape;
        } break;
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
    if(forced_collision_type!=-1)
        for(int i=0;i<collision_objects.m;i++)
            collision_objects(i)->type=(COLLISION_TYPE)forced_collision_type;
}
template class STANDARD_TESTS<VECTOR<float,2> >;
template class STANDARD_TESTS<VECTOR<double,2> >;
}
