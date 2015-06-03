//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Tools/Matrices/FRAME.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_OLDROYD_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/OLDROYD_NEO_HOOKEAN.h>
#include "STANDARD_TESTS_2D.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,2> >::
STANDARD_TESTS(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_BASE<TV>(stream_type,parse_args),Nsurface(0)
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
            VECTOR<T,1> angular_velocity(0.4*scale_speed);
            T density=2*scale_mass;
            Seed_Particles_Helper(sphere,[=](const TV& X){return angular_velocity.Cross(X-sphere.center);},
                [=](const TV&){return MATRIX<T,2>::Cross_Product_Matrix(angular_velocity);},
                density,particles_per_cell);
            T total_mass=particles.mass.Sum();
            TV total_momentum=particles.V.Weighted_Sum(particles.mass);
            TV dV=total_momentum/total_mass;
            particles.V-=dV;
            Add_Neo_Hookean(1e3*scale_E,0.3);
        } break;
        case 2:{ // oscillating circle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles_Helper(sphere,[=](const TV& X){return TV(0.1,0);},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>()+1.5);
            Add_Fixed_Corotated(1e3*scale_E,0.3);
        } break;
        case 3:{ // freefall circle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles_Helper(sphere,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            Add_Gravity(TV(0,-9.8));
        } break;
        case 4:{ // colliding of two rings
            if(!user_resolution) resolution=48;
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(),TV(0.48,0.48)),true);
            ARRAY<SPHERE<TV> > spheres; ARRAY<TV> v0; ARRAY<T> r;
            spheres.Append(SPHERE<TV>(TV(0.1,0.24),0.04)); spheres.Append(SPHERE<TV>(TV(0.4,0.24),0.04));
            v0.Append(TV(50,0)); v0.Append(TV(-50,0));
            r.Append(0.03); r.Append(0.03);
            for(int s=0;s<spheres.m;s++){
                SPHERE<TV>& sphere=spheres(s);
                T density=1010*scale_mass;
                int last=particles.number;
                Seed_Particles_Helper(sphere,[=](const TV& X){return v0(s);},[=](const TV&){return MATRIX<T,2>();},density,particles_per_cell);
                for(int k=last;k<particles.number;k++){
                    if((particles.X(k)-sphere.center).Magnitude_Squared()<sqr(r(s))){
                        particles.deletion_list.Append(k);}}}
            particles.Delete_Elements_On_Deletion_List();
            Add_Neo_Hookean(0.073e9*scale_E,0.4);
        } break;
        case 5:{ // rebound of an elastic cylinder
            if(!user_resolution) resolution=10;
            T dx=(T)5/resolution;
            grid.Initialize(TV_INT(3,1)*resolution+9,RANGE<TV>(TV(),TV(15,5)).Thickened(dx*(T)4.5),true);
            Add_Collision_Object(RANGE<TV>(TV(-5,-5),TV(0+dx/2,15)),COLLISION_TYPE::slip,0);
            Add_Collision_Object(RANGE<TV>(TV(15-dx/2,-5),TV(20,15)),COLLISION_TYPE::slip,0);
            SPHERE<TV> sphere(TV(2.5,2.5),1.5);
            T density=4*scale_mass;
            Seed_Particles_Helper(sphere,[=](const TV& X){return TV(0.5,0);},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            Add_Neo_Hookean(85.5*scale_E,0.425); //solve({E/(2*(1+r))=30,E*r/((1+r)*(1-2*r))=170},{E,r});
        } break;
        case 6:{ // skew impact of two elastic cylinders
            if(!user_resolution) resolution=12;
            int third_resolution=(resolution+2)/3;
            T dx=(T)4/third_resolution;
            grid.Initialize(TV_INT(5,3)*third_resolution+9,RANGE<TV>(TV(),TV(20,12)).Thickened(dx*(T)4.5),true);
            T density=5*scale_mass;
            SPHERE<TV> sphere1(TV(3,3),2);
            Seed_Particles_Helper(sphere1,[=](const TV& X){return TV(0.75,0);},[=](const TV&){return MATRIX<T,2>();},density,particles_per_cell);
            SPHERE<TV> sphere2(TV(16,5),2);
            Seed_Particles_Helper(sphere2,[=](const TV& X){return TV(-0.75,0);},[=](const TV&){return MATRIX<T,2>();},density,particles_per_cell);
            Add_Neo_Hookean(31.685*scale_E,0.44022); //solve({E/(2*(1+r))=11,E*r/((1+r)*(1-2*r))=81},{E,r});
        } break;
        case 7:{ // ping-pong ring
            // ./mpm 7 -flip 0  -affine -midpoint -max_dt 1e-3 -cfl .1 -framerate 2400 -newton_tolerance 1e-5 -solver_tolerance 1e-5  -last_frame 240 -order 2 -print_stats | grep 'total'
            if(!user_resolution) resolution=480;
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(),TV(0.48,0.48)),true);
            Add_Collision_Object(RANGE<TV>(TV(-5,-5),TV(0.11,15)),COLLISION_TYPE::separate,0);
            Add_Collision_Object(RANGE<TV>(TV(0.3,-5),TV(20,15)),COLLISION_TYPE::separate,0);
            ARRAY<SPHERE<TV> > spheres; ARRAY<TV> v0; ARRAY<T> r;
            spheres.Append(SPHERE<TV>(TV(0.2,0.24),0.04));
            v0.Append(TV(50,0));
            r.Append(0.03);
            for(int s=0;s<spheres.m;s++){
                SPHERE<TV>& sphere=spheres(s);
                T density=1010*scale_mass;
                int last=particles.number;
                Seed_Particles_Helper(sphere,[=](const TV& X){return v0(s);},[=](const TV&){return MATRIX<T,2>();},density,particles_per_cell);
                for(int k=last;k<particles.number;k++){
                    if((particles.X(k)-sphere.center).Magnitude_Squared()<sqr(r(s))){
                        particles.deletion_list.Append(k);}}}
            particles.Delete_Elements_On_Deletion_List();
            Add_Neo_Hookean(0.073e9*scale_E,0.4);
        } break;
        case 8:{ // collision an elastic cylinder (TODO: fix description.)
            if(!user_resolution) resolution=10;
            grid.Initialize(TV_INT()+resolution+9,RANGE<TV>(TV(),TV(5,5)),true);
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1,false);
            Add_Collision_Object(SPHERE<TV>(TV(4,3),1),COLLISION_TYPE::separate,.3);
            SPHERE<TV> sphere(TV(2.55,3.55),.3);
            T density=4*scale_mass;
            Seed_Particles_Helper(sphere,[=](const TV& X){return TV(3.0,0);},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            Add_Neo_Hookean(scale_E,0.425);
            Add_Gravity(TV(0,-1.8));
        } break;
        case 9:{ // collision an elastic cylinder (TODO: fix description.)
            if(!user_resolution) resolution=10;
            grid.Initialize(TV_INT()+resolution+9,RANGE<TV>(TV(),TV(5,5)),true);
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1,false);
            SPHERE<TV> sphere(TV(2.55,3.55),.3);
            T density=4*scale_mass;
            Seed_Particles_Helper(sphere,[=](const TV& X){return TV(3.0,0);},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            Add_Neo_Hookean(scale_E,0.425,&mpm_particles);

            SPHERE<TV> sphere2(TV(3.55,3.35),.3);
            TRIANGULATED_AREA<T>* ta=TESSELLATION::Generate_Triangles(sphere2,ceil(sphere2.radius/grid.dX.x*2));
            TRIANGULATED_AREA<T>& new_ta=Seed_Lagrangian_Particles(*ta,[=](const TV& X){return TV(-3.0,0);},
                [=](const TV&){return MATRIX<T,2>();},density,true);
            Add_Neo_Hookean(new_ta,scale_E,0.425);

            Add_Gravity(TV(0,-1.8));
        } break;
        case 10:{ // mpm projectile vs end-holded wall
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1,false);
            Add_Collision_Object(RANGE<TV>(TV(.45,.75),TV(.65,.85)),COLLISION_TYPE::stick,0);
            Add_Collision_Object(RANGE<TV>(TV(.45,.15),TV(.65,.25)),COLLISION_TYPE::stick,0);
            {SPHERE<TV> sphere(TV(.2,.5),.06);
                T density=2*scale_mass;
                Seed_Particles_Helper(sphere,[=](const TV& X){return TV(1,0);},[=](const TV&){return MATRIX<T,2>();},
                    density,particles_per_cell);
                ARRAY<int> foo(IDENTITY_ARRAY<>(particles.number));
                Add_Fixed_Corotated(20*scale_E,0.3,&foo);}
            int N_sphere_particles=particles.number;
            {RANGE<TV> box(TV(.5,.2),TV(.6,.8));
                T density=1*scale_mass;
                Seed_Particles_Helper(box,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},
                    density,particles_per_cell);
                int N_box_particles=particles.number-N_sphere_particles;
                ARRAY<int> foo(N_box_particles);
                for(int k=0;k<foo.m;k++) foo(k)=k+N_sphere_particles;
                Add_Fixed_Corotated(1*scale_E,0.3,&foo);}

            Add_Gravity(TV(0,-1.8));
        } break;
        case 11:{ // freefall circle, rising ground
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles_Helper(sphere,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            ARRAY<int> foo(IDENTITY_ARRAY<>(particles.number));
            Add_Fixed_Corotated(1*scale_E,0.3,&foo);
            Add_Gravity(TV(0,-9.8));
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1,false);
            MPM_COLLISION_IMPLICIT_OBJECT<TV>* bottom=dynamic_cast<MPM_COLLISION_IMPLICIT_OBJECT<TV>*>(collision_objects(3));
            bottom->func_frame=[this](T time){return FRAME<TV>(TV(0,(T).75-abs((T).15*time*scale_speed-(T).75)));};
            bottom->func_twist=[this](T time){return TWIST<TV>(TV(0,-sign((T).15*time*scale_speed-(T).75)*(T).15*scale_speed),typename TV::SPIN());};
        } break;
        case 12:{ // freefall circle, rising ground, penalty collisions
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles_Helper(sphere,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            ARRAY<int> foo(IDENTITY_ARRAY<>(particles.number));
            Add_Fixed_Corotated(1*scale_E,0.3,&foo);
            Add_Gravity(TV(0,-9.8));
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1,true);
        } break;
        case 13:{ // surface tension test: fixed topology circle shell
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SEGMENTED_CURVE_2D<T>* sc=SEGMENTED_CURVE_2D<T>::Create();
            TV center(0.5,0.5);
            T radius=0.2;
            int N=80;
            T density=1*scale_mass;
            sc->particles.Add_Elements(N);
            sc->Update_Number_Nodes();
            for(int n=0;n<N;n++){
                T theta=(T)n*2.0*3.1415926/(T)N;
                TV X=center+TV(radius*cos(theta),radius*sin(theta));
                X(0)+=random.Get_Uniform_Number(-radius*0.2,radius*0.2);
                X(1)+=random.Get_Uniform_Number(-radius*0.2,radius*0.2);
                LOG::cout<<X<<std::endl;
                sc->particles.X(n)=X;
                int np1=(n!=N-1)?(n+1):0;
                sc->mesh.elements.Append(TV_INT(n,np1));
                LOG::cout<<sc->mesh.elements(n)<<std::endl;}
            SEGMENTED_CURVE_2D<T>& new_sc=Seed_Lagrangian_Particles(*sc,[=](const TV& X){return TV(0.0,0);},[=](const TV&){return MATRIX<T,2>();},density,true);
            SURFACE_TENSION_FORCE<TV>* stf=new SURFACE_TENSION_FORCE<TV>(new_sc,(T)0.01);
            stf->use_velocity_independent_implicit_forces=true;
            Add_Force(*stf);
        } break;
        case 14:{ // test dynamic changing lagrangian mesh (in Begin_Frame)
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(-1,-1),TV(2,2)),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles_Helper(sphere,[=](const TV& X){return TV(0,0);},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>()+1.5);
            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            Add_Fixed_Corotated(scale_E,0.3,&mpm_particles);
        } break;
        case 15:{ // surface tension circle

            //one   ./mpm 15 -resolution 64 -last_frame 500 -scale_E 10 -newton_tolerance 1e-5 -particles_per_cell 8 -max_dt 1e-3 -affine

            //two ./mpm 15   -resolution 64 -last_frame 500 -scale_E 10 -newton_tolerance 1e-5 -particles_per_cell 4  -max_dt 1e-3 -affine

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            T density=2*scale_mass;
            if(1) // c=1e-2
            {
                RANGE<TV> box(TV(.3,.3),TV(.5,.5));
                Seed_Particles_Helper(box,[=](const TV& X){return TV(0,0);},[=](const TV&){return MATRIX<T,2>();},
                    density,particles_per_cell);
                RANGE<TV> box2(TV(.5,.35),TV(.8,.45));
                Seed_Particles_Helper(box2,[=](const TV& X){return TV(0,0);},[=](const TV&){return MATRIX<T,2>();},
                    density,particles_per_cell);
                ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
                bool no_mu=true;
                Add_Fixed_Corotated(scale_E,0.3,&mpm_particles,no_mu);
            }
            if(0) // c=1e-3
            {
                RANGE<TV> box(TV(.2,.3),TV(.4,.5));
                Seed_Particles_Helper(box,[=](const TV& X){return TV(0.1,0);},[=](const TV&){return MATRIX<T,2>();},
                    density,particles_per_cell);
                RANGE<TV> box2(TV(.6,.3),TV(.8,.5));
                Seed_Particles_Helper(box2,[=](const TV& X){return TV(-0.1,0);},[=](const TV&){return MATRIX<T,2>();},
                    density,particles_per_cell);
                ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
                bool no_mu=true;
                Add_Fixed_Corotated(scale_E,0.3,&mpm_particles,no_mu);
            }
        } break;
        case 16:{ // oscillating circle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            use_oldroyd=true;
            Seed_Particles_Helper(sphere,[=](const TV& X){return TV(0.1,0);},[=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>()+1.5);
            particles.S.Fill(SYMMETRIC_MATRIX<T,2>()+sqr(1.5));

            OLDROYD_NEO_HOOKEAN<TV> *neo=new OLDROYD_NEO_HOOKEAN<TV>;
            neo->mu=10000;
            neo->lambda=10000;
            MPM_OLDROYD_FINITE_ELEMENTS<TV> *fe=new MPM_OLDROYD_FINITE_ELEMENTS<TV>(force_helper,*neo,gather_scatter,0);
            Add_Force(*fe);
        } break;
        case 21:{ // circle with random initial velocities
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(-3,-3),TV(4,4)),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles_Helper(sphere,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},density,particles_per_cell);
            random.Fill_Uniform(particles.V,-1,1);
            Add_Fixed_Corotated(1e3*scale_E,0.3);
        } break;
        case 22:{ // (fluid test) pool of water 
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(grid.dX*(T).5,TV(1-grid.dX(0)*(T).5,0.25));
            T density=2*scale_mass;
            Seed_Particles_Helper(box,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},density,particles_per_cell);
            Add_Fixed_Corotated(1e3*scale_E,0.3);
            Add_Gravity(TV(0,-1.8));
        } break;
        case 23:{ // (fluid test) dam break 
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box;
            box=RANGE<TV>(grid.dX*(T).5,TV(0.2,0.75));
            T density=2*scale_mass;
            Seed_Particles_Helper(box,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},density,particles_per_cell);
            Add_Fixed_Corotated(1e3*scale_E,0.3);
            Add_Gravity(TV(0,-1.8));
        } break;
        case 24:{ // (fluid test) circle drop 
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.7),.2);
            T density=2*scale_mass;
            Seed_Particles_Helper(sphere,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},density,particles_per_cell);
            Add_Fixed_Corotated(1e3*scale_E,0.3);
            Add_Gravity(TV(0,-1.8));
        } break;
        case 25:{ // (fluid test) pool of water w/ single particle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(grid.dX*(T).5,TV(1-grid.dX(0)*(T).5,0.25));
            T density=2*scale_mass;
            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            Seed_Particles_Helper(box,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},density,particles_per_cell);
            Add_Particle(TV(.5,.9),TV(),mass,volume,MATRIX<T,TV::m>()+1,MATRIX<T,TV::m>());
            Add_Gravity(TV(0,-1.8));
        } break;
        case 26:{ // Raleigh Taylor
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(grid.dX*(T).5,TV(1-grid.dX(0)*(T).5,0.20));
            T density=2*scale_mass;
            Seed_Particles_Helper(box,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},density,particles_per_cell);
            density*=10;
            box+=TV(0,0.20);
            Seed_Particles_Helper(box,[=](const TV& X){return TV();},[=](const TV&){return MATRIX<T,2>();},density,particles_per_cell);
            Add_Gravity(TV(0,-1.8));
        } break;
        case 99:{ // single particle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            T density=2*scale_mass;
            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            Add_Particle(TV(.8,.5),TV(),mass,volume,MATRIX<T,TV::m>()+1,MATRIX<T,TV::m>());
            Add_Gravity(TV(0,-1.8));
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
    switch(test_number)
    {
        case 14:{
            if(frame==0){
                T density=2*scale_mass;
                SPHERE<TV> sphere2(TV(0.55,0.35),.2);
                TRIANGULATED_AREA<T>* ta=TESSELLATION::Generate_Triangles(sphere2,ceil(sphere2.radius/grid.dX.x*2));
                Nsurface=ta->particles.number;
                TRIANGULATED_AREA<T>& new_ta=Seed_Lagrangian_Particles(*ta,[=](const TV& X){return TV(-3.0,0);},
                    [=](const TV&){return MATRIX<T,2>();},density,true);
                Add_Fixed_Corotated(new_ta,(T)10*scale_E,0.425);}
            else if(frame==6){
                for(int k=particles.number-Nsurface;k<particles.number;k++) particles.Add_To_Deletion_List(k);
                particles.Delete_Elements_On_Deletion_List();
                T density=2*scale_mass;
                SPHERE<TV> sphere2(TV(1,1),.3);
                TRIANGULATED_AREA<T>* ta=TESSELLATION::Generate_Triangles(sphere2,ceil(sphere2.radius/grid.dX.x*2));
                Nsurface=ta->particles.number;
                TRIANGULATED_AREA<T>& new_ta=Seed_Lagrangian_Particles(*ta,[=](const TV& X){return TV(-3.0,-1.5);},
                    [=](const TV&){return MATRIX<T,2>();},density,true);
                Add_Fixed_Corotated(new_ta,(T)10*scale_E,0.425);}
        } break;
    }
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
    // static int count=0;
    switch(test_number)
    {
        case 15:{
            // if(count++>1 ) break;
            int N_non_surface=particles.number-Nsurface;
            for(int k=N_non_surface;k<particles.number;k++){
                particles.Add_To_Deletion_List(k);

                // return stolen property
                int m=steal(k-N_non_surface);
                particles.mass(m)+=particles.mass(k);
                particles.volume(m)+=particles.volume(k);
                // particles.X(m)=(particles.X(m)+particles.X(k))*0.5;
                particles.V(m)=(particles.V(m)+particles.V(k))*0.5;
                if(use_affine) particles.B(m)=(particles.B(m)+particles.B(k))*0.5;


            }
            LOG::cout<<"deleting "<<Nsurface<<" particles..."<<std::endl;
            particles.Delete_Elements_On_Deletion_List();
            lagrangian_forces.Delete_Pointers_And_Clean_Memory();
            this->deformable_body_collection.structures.Delete_Pointers_And_Clean_Memory();

            SEGMENTED_CURVE_2D<T>* surface=SEGMENTED_CURVE_2D<T>::Create();
            MARCHING_CUBES<TV>::Create_Surface(*surface,grid,mass,particles.mass(0)*1);

            // blobby
            if(0){
                T r=grid.DX()(0)*.9;
                GRID<TV> lsgrid; lsgrid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box());

                ARRAY<T,TV_INT> phi;
                phi.Resize(lsgrid.Domain_Indices(ghost));
                phi.Fill(FLT_MAX);

                for(int k=0;k<particles.number;k++){
                    TV Xp=particles.X(k);
                    TV_INT cell=lsgrid.Cell(Xp,ghost);
                    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT()-3,TV_INT()+7));it.Valid();it.Next()){
                        TV_INT index=cell+it.index;
                        TV Xi=lsgrid.Node(index);
                        T d=(Xi-Xp).Magnitude();
                        if(d<phi(index)) phi(index)=d;}}
                for(int i=0;i<phi.array.m;i++) phi.array(i)-=r;

                MARCHING_CUBES<TV>::Create_Surface(*surface,lsgrid,phi,0);
            } // end blobby

            T min_marching_length=FLT_MAX;
            for(int k=0;k<surface->mesh.elements.m;k++){
                TV A=surface->particles.X(surface->mesh.elements(k).x);
                TV B=surface->particles.X(surface->mesh.elements(k).y);
                T d=(A-B).Magnitude();
                min_marching_length=min(min_marching_length,d);}
            LOG::cout<<"MIN MARCHING LENGTH: " <<min_marching_length<<std::endl;

            int Nold=particles.number;
            Nsurface=surface->particles.number;
            LOG::cout<<"adding "<<Nsurface<<" particles..."<<std::endl;
            SEGMENTED_CURVE_2D<T>& new_sc=Seed_Lagrangian_Particles(*surface,[=](const TV& X){return TV(0.0,0);},[=](const TV&){return MATRIX<T,2>();},(T)0.001,true);

            // for(int k=Nold;k<particles.number;k++) particles.mass(k)=surface_mass/Nsurface;

            steal.Clean_Memory();
            for(int k=Nold;k<particles.number;k++){
                
                T dist=FLT_MAX;
                int m=-1;
                for(int q=0;q<Nold;q++){
                    if(steal.Contains(q)) continue;
                    T dd=(particles.X(q)-particles.X(k)).Magnitude_Squared();
                    if(dd<dist){
                        dist=dd;
                        m=q;}}
                T split_mass=particles.mass(m)*.5;
                T split_volume=particles.volume(m)*.5;
                TV com=particles.X(m);

                particles.mass(m)=split_mass;
                particles.volume(m)=split_volume;
                // particles.X(m)=com*(T)2-particles.X(k);
                
                particles.mass(k)=split_mass;
                particles.volume(k)=split_volume;
                particles.V(k)=particles.V(m);
                if(use_affine) particles.B(k)=particles.B(m);

                steal.Append(m);

            }


            SURFACE_TENSION_FORCE<TV>* stf=new SURFACE_TENSION_FORCE<TV>(new_sc,(T)1e-2);
            stf->use_velocity_independent_implicit_forces=true;
            Add_Force(*stf);

            // Dump_Surface(new_sc,VECTOR<T,3>(1,1,0)); 
        } break;
    }

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
