//#####################################################################
// Copyright 2015, Lin Huang, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INTERSECTION.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INVERT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_UNION.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include "STANDARD_TESTS_2D.h"
#include <fstream>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,2> >::
STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_BASE<TV>(stream_type_input,parse_args)
{
    parse_args.Parse();
    if(!this->override_output_directory) viewer_dir.output_directory=LOG::sprintf("Test_%i",test_number);
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
        side_bc_type.Fill(BC_PERIODIC);
    switch(test_number)
    {
        case 1:{ // stationary circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            density=unit_rho*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
        } break;
        case 2:
        case 15:{ // translating circle
            if(test_number==2) side_bc_type.Fill(BC_FREE);
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            density=unit_rho*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV(m/s,0);},0,density,particles_per_cell);
        } break;
        case 3:{ // rotating circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            VECTOR<T,1> angular_velocity(0.4/s);
            density=unit_rho*scale_mass;
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
            density=unit_rho*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
        } break;
        case 5:{ // stationary pool
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            density=unit_rho*scale_mass;
            Seed_Particles(RANGE<TV>(TV(),TV(m,.5*m)),0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
        } break;
        case 6:{ // freefall rectangle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            density=unit_rho*scale_mass;
            Seed_Particles(RANGE<TV>(TV(.2*m,.2*m),TV(.5*m,.8*m)),0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
        } break;
        case 100: this->Commandline_Analytic_test(); break;
        case 8:{ // concave shape
            this->use_phi=true;
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            density=unit_rho*scale_mass;
            SPHERE<TV> sphere0(TV(.3,.3)*m,.1*m);
            SPHERE<TV> sphere1(TV(.44,.44)*m,.1*m);
            auto shape=Unite(Make_IO(sphere0),Make_IO(sphere1));
            Seed_Particles(*shape,0,0,density,particles_per_cell);
            delete shape;
        } break;
        case 11:
        case 12:{ // free fall circle with curved boundary
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            if(test_number==12) sphere.radius=.4*m;
            density=unit_rho*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
            // add a circle collision object
            SPHERE<TV> sph(TV(.5,.5)*m,.4*m);
            Add_Collision_Object(Invert(Make_IO(sph)),COLLISION_TYPE::slip,0,0,0);
        } break;
        case 13:{ // half filled stationary pool with curved boundary
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            density=unit_rho*scale_mass;
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
            side_bc_type.Fill(BC_PERIODIC);
            VECTOR<T,1> angular_velocity(0.4/s);
            density=unit_rho*scale_mass;
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
            side_bc_type.Fill(BC_PERIODIC);
            T a=extra_int.m>=1?extra_int(0):1;
            Set_Grid(RANGE<TV>::Centered_Box()*pi*m);
            density=unit_rho*scale_mass;
            T b=-2*viscosity/density*sqr(a)/s;
            a/=m;
            use_analytic_field=true;
            Add_Velocity([=](auto X,auto t){return (rot*cos(X*a))*sin(X*a)*exp(b*t)*m/s;});
            Add_Pressure([=](auto X,auto t){return this->density/4*exp(2*b*t)*(cos(2*a*X(0))+cos(2*a*X(1)))*sqr(m/s);});
            Setup_Analytic_Boundary_Conditions();
            Seed_Particles_Analytic(grid.domain,density,particles_per_cell);
            end_frame.Append([=](int frame){Check_Analytic_Velocity();});
        } break;
        case 18:{ // full of fluid with random inital velocities
            T a=extra_T.m>=1?extra_T(0):-1;
            T b=extra_T.m>=2?extra_T(1):1;
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            density=unit_rho*scale_mass;
            RANDOM_NUMBERS<T> rand;
            rand.Set_Seed();
            auto V_func=[&](const TV&){
                return TV(rand.Get_Uniform_Number(a,b),rand.Get_Uniform_Number(a,b))*m/s;};
            Seed_Particles(grid.domain,V_func,0,density,particles_per_cell);
        } break;
        case 21:{ // vortex shedding test
            density=unit_rho*scale_mass;
            T velocity=1;
            Set_Grid(RANGE<TV>(TV(-2,-4),TV(14,4))*m,TV_INT(2,1));
            SPHERE<TV> sphere(TV(),.25*m);
            auto shape=Intersect(Make_IO(grid.domain),Invert(Make_IO(sphere)));
            Seed_Particles(*shape,0,0,density,particles_per_cell);
            delete shape;
            side_bc_type(1)=BC_FREE;
            bc_velocity=[=](const TV& X,T time){return TV(X.x<-1.999?velocity:0,0);};
            Add_Collision_Object(sphere,COLLISION_TYPE::slip,0,0,0);
            RANGE<TV> source_range=grid.domain;
            source_range.min_corner.x-=grid.dX.x;
            Add_Source(grid.domain.min_corner,TV(1,0),TV(velocity,0),TV(),
                Make_IO(source_range),density,particles_per_cell,true);
            auto write=[this]()
            {
                static int id=0;
                Create_Directory(viewer_dir.output_directory+LOG::sprintf("/v%d",id));
                Write_To_File(stream_type,LOG::sprintf("%s/v%d/mac_velocities",viewer_dir.output_directory.c_str(),id),this->velocity);
                Write_To_File(stream_type,LOG::sprintf("%s/v%d/mass",viewer_dir.output_directory.c_str(),id),mass);
                Write_To_File(stream_type,LOG::sprintf("%s/v%d/psi_N",viewer_dir.output_directory.c_str(),id),this->psi_N);
                Write_To_Text_File(LOG::sprintf("%s/v%d/time",viewer_dir.output_directory.c_str(),id),time,"\n");
                Write_To_Text_File(viewer_dir.output_directory+"/common/last_grid_data",id,"\n");
                id++;
            };
            Add_Callbacks(false,"p2g",[write](){
                    static bool done=false;
                    if(!done) write();
                    done=true;});
            Add_Callbacks(false,"time-step",write);
        } break;
        case 22:{ // Analytic velocity field generated by streamfunction xy(1-x)(1-y)(xy+yt^2-2xt)
            side_bc_type.Fill(BC_SLIP);
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            density=unit_rho*scale_mass;
            use_analytic_field=true;
            auto f=[=](auto x){return x*(1-x)*(x*x-x-1);};
            auto df=[=](auto x){return (1-2*x)*(2*x*x-2*x-1);};
            auto g=[=](auto x){return x*(1-x)*(x+1)*(3*x*x-7);};
            auto dg=[=](auto x){return (1-3*x*x)*(3*x*x-7)+x*(1-x)*(x+1)*(6*x);};
            int a=1,b=0,c=1;
            Add_Velocity(
                [=](auto Y,auto t)
                {
                    auto X=Y/m;
                    auto x=X(0),y=X(1);
                    return (T)a*Make_Vector<T>(-df(y)*f(x),df(x)*f(y))+(b+c*t)*Make_Vector<T>(-df(y)*g(x),dg(x)*f(y));
                });
            Add_Pressure([](auto X,auto t){
                    auto x=X(0),y=X(1);
                    auto p=x*y*(1-x)*(1-y);
                    return (X(0)-X(0)*X(1)+X(1)*X(1)+t)*p;
                });
            Setup_Analytic_Boundary_Conditions();
            Seed_Particles_Analytic(grid.domain,density,particles_per_cell);
            end_frame.Append([=](int frame){Check_Analytic_Velocity();});
        } break;
        case 23:{ // Dirichlet BC convergence test
            Set_Grid(RANGE<TV>::Unit_Box());
            density=unit_rho*scale_mass;
            use_analytic_field=true;
            Add_Velocity([=](auto X,auto t){return Auto_Hess_Vector(X(0),-X(1))*m/s;});
            Add_Pressure([=](auto X,auto t){return (X(0)-X(0)*X(1)+X(1)*X(1)+t)*unit_p;});
            Setup_Analytic_Boundary_Conditions();
            SPHERE<TV> sphere(TV(.5,.5),.1);
            Seed_Particles_Analytic(sphere,density,particles_per_cell);
            auto valid_face=[=](const FACE_INDEX<TV::m>& face)
            {
                TV X=grid.Face(face);
                T cx=X(0)/exp(time),cy=X(1)/exp(-time);
                return sqr(cx-.5)+sqr(cy-.5)<=sqr(.1);
            };
            end_frame.Append([=](int frame){Check_Analytic_Velocity(valid_face);});
        } break;
        case 24:{ // Analytic velocity field generated by streamfunction (xy(1-x)(1-y))^n
            //side_bc_type.Fill(BC_NOSLIP);
            side_bc_type.Fill(BC_PERIODIC);
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            density=unit_rho*scale_mass;
            use_analytic_field=true;
            Add_Velocity(
                [=](auto Y,auto t)
                {
                    auto X=Y/m;
                    TV o(1,1);
                    auto x=X(0),y=X(1);
                    auto d=(1-x)*x*(1-y)*y;
                    auto dd=Make_Vector<T>((1-2*x)*(1-y)*y,(1-x)*x*(1-2*y));
                    auto p=x+2*y+1;
                    TV dp(1,2);
                    return d*d*rot*(d*dp+3*p*dd);
                });
            Add_Pressure([](auto X,auto t){
                    auto x=X(0),y=X(1);
                    auto p=x*y*(1-x)*(1-y);
                    return (X(0)-X(0)*X(1)+X(1)*X(1)+t)*p;
                });
            Setup_Analytic_Boundary_Conditions();
            Seed_Particles_Analytic(grid.domain,density,particles_per_cell);
            end_frame.Append([=](int frame){Check_Analytic_Velocity();});
        } break;
        case 25:{ // Dam break
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            density=unit_rho*scale_mass;
            Seed_Particles(RANGE<TV>(TV(.0,.0),TV(m/3,2*m/3)),0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
        } break;
        case 26:{ // One Taylor-Green Vortex in [0,pi]^2, with slip BC.
            side_bc_type.Fill(BC_SLIP);
            Set_Grid(RANGE<TV>::Unit_Box()*pi*m);
            density=unit_rho*scale_mass;
            T b=-2*viscosity/density/s;
            use_analytic_field=true;
            Add_Velocity([=](auto X,auto t){return (rot*cos(X/m))*sin(X/m)*exp(b*t)*m/s;});
            Add_Pressure([=](auto X,auto t){return X(0)-X(0)*X(1)+X(1)*X(1)+t;});
            Setup_Analytic_Boundary_Conditions();
            Seed_Particles_Analytic(grid.domain,density,particles_per_cell);
            end_frame.Append([=](int frame){Check_Analytic_Velocity();});
        } break;
        case 27:{ // Periodic version of "square" test (22)
            side_bc_type.Fill(BC_PERIODIC);
            Set_Grid(RANGE<TV>::Centered_Box()*m);
            density=unit_rho*scale_mass;
            use_analytic_field=true;
            Add_Velocity(
                [=](auto Y,auto t)
                {
                    auto X=Y/m;
                    TV o(1,1);
                    return X*(o-abs(X))*(rot*(o-abs(X)*2));
                });
            Add_Pressure([](auto X,auto t){
                    //auto x=X(0),y=X(1);
                    //auto p=x*y*(1-x)*(1-y);
                    //return (X(0)-X(0)*X(1)+X(1)*X(1)+t)*p;
                    return 0;
                });
            Setup_Analytic_Boundary_Conditions();
            Seed_Particles_Analytic(grid.domain,density,particles_per_cell);
            end_frame.Append([=](int frame){Check_Analytic_Velocity();});
        } break;
        case 28:{ // inlet-cube
            density=unit_rho*scale_mass;
            T velocity=extra_T.m>=1?extra_T(0):0.25,stop_time=extra_T.m>=2?extra_T(1):200*frame_dt;
            Set_Grid(RANGE<TV>::Unit_Box());
            Seed_Particles(grid.domain,0,0,density,particles_per_cell);
            T inlet_center=0.625,inlet_radius=0.125;
            this->bc_type=[=](const TV& X,T time)
            {
                if(X.y>0.001*grid.dX(1) || time>stop_time) return BC_SLIP;
                else if(abs(X.x-0.2)<=0.125/2) return BC_FREE;
                else if(abs(X.x-0.9)<=0.125/2) return BC_FREE;
                else return BC_SLIP;
            };
            bc_velocity=[=](const TV& X,T time)
            {
                return TV(0,(time<=stop_time&&abs(X.x-inlet_center)<=inlet_radius&&X.y<0.001*grid.dX(1))?velocity:0);
            };
            RANGE<TV> source_range=grid.domain;
            source_range.min_corner.y-=.5*grid.dX.y;
            source_range.min_corner.x=inlet_center-inlet_radius;
            source_range.max_corner.x=inlet_center+inlet_radius;
            source_range.max_corner.y=0;
            TV X0=source_range.min_corner;
            X0.x+=2*inlet_radius;
            Add_Source(X0,TV(0,1),TV(0,velocity),TV(),
                Make_IO(source_range),density,particles_per_cell,true);
        } break;
        case 29:{
            Set_Grid(RANGE<TV>::Centered_Box()*2);
            side_bc_type.Fill(BC_SLIP);
            //side_bc_type(1)=BC_FREE;
            density=unit_rho*scale_mass;
            use_analytic_field=true;
            Add_Velocity([=](auto X,auto t)
            {
                auto x=X(0),y=X(1);
                auto x2=x*x;
                auto y2=y*y;
                auto x4=x2*x2;
                auto y4=y2*y2;
                auto x6=x2*x4;
                auto y6=y2*y4;
                auto x8=x4*x4;
                auto y8=y4*y4;
                return Make_Vector<T>(
                    -(T)2/1775*(x-2)*(x+2)*y*
                    (5325*x8*y4+14200*x6*y6+8875*x4*y8-43370*x8*y2-185367*x6*y4-218756*x4*y6
                     -72925*x2*y8+94440*x8+799054*x6*y2+1615671*x4*y4+1070044*x2*y6+180500*y8
                     -1120888*x6-4676038*x4*y2-5345697*x2*y4-1714416*y6+4413464*x4+10176970*x2*y2
                     +5308260*y4-5828392*x2-6256616*y2+2604960),
                    (T)2/1775*(y-2)*(y+2)*x*
                    (8875*x8*y4+14200*x6*y6+5325*x4*y8-72925*x8*y2-218756*x6*y4-185367*x4*y6
                     -43370*x2*y8+180500*x8+1070044*x6*y2+1615671*x4*y4+799054*x2*y6+94440*y8
                     -1714416*x6-5345697*x4*y2-4676038*x2*y4-1120888*y6+5308260*x4+10176970*x2*y2
                     +4413464*y4-6256616*x2-5828392*y2+2604960))/(T)1e5;
            });
            Add_Pressure([=](auto X,auto t)
            {
                auto x=X(0),y=X(1);
                return x-x*y+y*y+t;
            });
            Setup_Analytic_Boundary_Conditions();
            SPHERE<TV> sphere(TV(),1);
            auto shape=Intersect(Make_IO(grid.domain),Invert(Make_IO(sphere)));
            Seed_Particles_Analytic(*shape,density,particles_per_cell);
            //Seed_Particles_Analytic(grid.domain,density,particles_per_cell);
            delete shape;
            Add_Collision_Object(sphere,COLLISION_TYPE::slip,0,0,0);
            end_frame.Append([=](int frame){Check_Analytic_Velocity();});
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
