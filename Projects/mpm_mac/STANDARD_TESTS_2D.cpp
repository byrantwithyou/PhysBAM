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
#include <Hybrid_Methods/Seeding/MPM_PARTICLE_SOURCE.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
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
    T m=this->m,s=this->s,unit_p=this->unit_p;
    MATRIX<T,2> rot(0,1,-1,0);
    if(bc_periodic)
        bc_type.Fill(BC_PERIODIC);
    switch(test_number)
    {
        case 1:{ // stationary circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=unit_rho*scale_mass;
            Set_Phases({density});
            Seed_Particles(sphere,0,0,density,particles_per_cell);
        } break;
        case 2:
        case 15:{ // translating circle
            if(test_number==2) bc_type.Fill(BC_INVALID);
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=unit_rho*scale_mass;
            Set_Phases({density});
            Seed_Particles(sphere,[=](const TV& X){return TV(m/s,0);},0,density,particles_per_cell);
        } break;
        case 3:{ // rotating circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            VECTOR<T,1> angular_velocity(0.4/s);
            T density=unit_rho*scale_mass;
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
            T density=unit_rho*scale_mass;
            Set_Phases({density});
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
        } break;
        case 5:{ // stationary pool
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=unit_rho*scale_mass;
            Set_Phases({density});
            Seed_Particles(RANGE<TV>(TV(),TV(m,.5*m)),0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
        } break;
        case 6:{ // freefall rectangle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=unit_rho*scale_mass;
            Set_Phases({density});
            Seed_Particles(RANGE<TV>(TV(.2*m,.2*m),TV(.5*m,.8*m)),0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
        } break;
        case 7:{ // stationary circles in two phases
            T density=unit_rho*scale_mass;
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
            T density=unit_rho*scale_mass;
            Set_Phases({density});
            SPHERE<TV> sphere0(TV(.3,.3)*m,.1*m);
            SPHERE<TV> sphere1(TV(.44,.44)*m,.1*m);
            auto shape=Unite(Make_IO(sphere0),Make_IO(sphere1));
            Seed_Particles(*shape,0,0,density,particles_per_cell);
            delete shape;
        } break;
        case 9: // freefall circles with different phases
        case 10:{ // freefall circles with same phase
            T density=unit_rho*scale_mass;
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
            T density=unit_rho*scale_mass;
            Set_Phases({density});
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
            // add a circle collision object
            SPHERE<TV> sph(TV(.5,.5)*m,.4*m);
            Add_Collision_Object(Invert(Make_IO(sph)),COLLISION_TYPE::slip,0,0,0);
        } break;
        case 13:{ // half filled stationary pool with curved boundary
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=unit_rho*scale_mass;
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
            bc_type.Fill(BC_PERIODIC);
            VECTOR<T,1> angular_velocity(0.4/s);
            T density=unit_rho*scale_mass;
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
        } break;
        case 16:{ // Initialize velocity field as Taylor-Green Vortex
            // To specify parameters, use: -I a
            bc_type.Fill(BC_PERIODIC);
            T a=extra_int.m>=1?extra_int(0):1;
            Set_Grid(RANGE<TV>::Centered_Box()*pi*m);
            T density=unit_rho*scale_mass;
            Set_Phases({density});
            T b=-2*mu/density*sqr(a)/s;
            a/=m;
            use_analytic_field=true;
            Add_Velocity([=](auto X,auto t){return (rot*cos(X*a))*sin(X*a)*exp(b*t)*m/s;});
            Add_Pressure([=](auto X,auto t){return density/4*exp(2*b*t)*(cos(2*a*X(0))+cos(2*a*X(1)))*sqr(m/s);});
            Setup_Analytic_Boundary_Conditions();
            Seed_Particles_Analytic(grid.domain,PHASE_ID(0),density,particles_per_cell);
            end_frame.Append([=](int frame){Check_Analytic_Velocity();});
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
            T density=unit_rho*scale_mass;
            Set_Phases({density});
            RANDOM_NUMBERS<T> rand;
            rand.Set_Seed();
            auto V_func=[&](const TV&){
                return TV(rand.Get_Uniform_Number(a,b),rand.Get_Uniform_Number(a,b))*m/s;};
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
            delete shape;
            delete cshape;
        } break;
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
        case 21:{ // vortex shedding test
            T density=unit_rho*scale_mass;
            T velocity=1;
            Set_Phases({density});
            particles.Store_Phase(true);
            Set_Grid(RANGE<TV>(TV(-2,-2),TV(14,2))*m,TV_INT(4,1));
            SPHERE<TV> sphere(TV(),m);
            auto shape=Intersect(Make_IO(grid.domain),Invert(Make_IO(sphere)));
            Seed_Particles(*shape,0,0,density,particles_per_cell);
            delete shape;
            bc_type(1)=BC_INVALID;
            bc_velocity(0)=[=](const TV& X,int axis,PHASE_ID pid,T time){return axis==0?velocity:0;};
            Add_Collision_Object(sphere,COLLISION_TYPE::slip,0,0,0);
            RANGE<TV> source_range=grid.domain;
            source_range.min_corner.x-=grid.dX.x;
            Add_Source(grid.domain.min_corner,TV(1,0),Make_IO(source_range),
                [=](TV X,T ts,T t,SOURCE_PATH<TV>& p)
                {
                    p.vp=TV(velocity,0);
                    p.xp=X+(t-ts)*p.vp;
                    p.xp_ts=-p.vp;
                    p.vp_ts=TV();
                    p.xp_x=MATRIX<T,TV::m>()+1;
                    p.vp_x=MATRIX<T,TV::m>();
                },density,particles_per_cell,true);
        } break;
        case 22:{ // Analytic velocity field generated by streamfunction xy(1-x)(1-y)(xy+yt^2-2xt)
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=unit_rho*scale_mass;
            Set_Phases({density});
            use_analytic_field=true;
            typedef DIAGONAL_MATRIX<T,2> TM;
            Add_Velocity(
                [=](auto Y,auto t)
                {
                    auto X=Y/m;
                    TV o(1,1);
                    auto x=X(0),y=X(1);
                    auto p=x*y*(1-x)*(1-y);
                    auto q=x*y+y*t*t-x*t*2;
                    auto dp=X*(X-o)*(X*2-o);
                    auto dq=Auto_Hess_Vector(-2*t+y,t*t+x);
                    auto dr=p*dq+q*dp;
                    return rot*dr;
                });
            Add_Pressure([](auto X,auto t){return X(0)-X(0)*X(1)+X(1)*X(1)+t;});
            Setup_Analytic_Boundary_Conditions();
            Seed_Particles_Analytic(grid.domain,PHASE_ID(0),density,particles_per_cell);
            end_frame.Append([=](int frame){Check_Analytic_Velocity();});
        } break;
        case 23:{ // Analytic velocity field generated by streamfunction xy+yt^2-2xt, with Dirichlet BC.
            Set_Grid(RANGE<TV>::Unit_Box());
            T density=unit_rho*scale_mass;
            Set_Phases({density});
            use_analytic_field=true;
            Add_Velocity([=](auto X,auto t){return Auto_Hess_Vector(-t*t/(s*s)-X(0)/m,-(T)2*t/s+X(1)/m)*m/s;});
            Add_Pressure([=](auto X,auto t){return (X(0)/s-X(0)*X(1)/m/s+X(1)*X(1)/m/s+t*m)*unit_p;});
            Setup_Analytic_Boundary_Conditions();
            SPHERE<TV> sphere(TV(.5,.5),.1);
            Seed_Particles_Analytic(sphere,PHASE_ID(0),density,particles_per_cell);
            end_frame.Append([=](int frame){Check_Analytic_Velocity();});
        } break;
    }
    if(mu){
        phases(PHASE_ID()).viscosity=mu;
        use_viscosity=true;}
    if(forced_collision_type!=-1)
        for(int i=0;i<collision_objects.m;i++)
            collision_objects(i)->type=(COLLISION_TYPE)forced_collision_type;
}
template class STANDARD_TESTS<VECTOR<float,2> >;
template class STANDARD_TESTS<VECTOR<double,2> >;
}
