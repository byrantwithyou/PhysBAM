//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Data_Structures/KD_TREE.h>
#include <Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Tools/Matrices/FRAME.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/PINNING_FORCE.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_OLDROYD_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_VISCOSITY.h>
#include <Hybrid_Methods/Forces/OLDROYD_NEO_HOOKEAN.h>
#include <Hybrid_Methods/Forces/VOLUME_PRESERVING_OB_NEO_HOOKEAN.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <fstream>
#include "STANDARD_TESTS_2D.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,2> >::
STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_BASE<TV>(stream_type_input,parse_args),
    use_surface_tension(false),Nsurface(0)
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
            Seed_Particles(sphere,[=](const TV& X){return angular_velocity.Cross(X-sphere.center);},
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
            Seed_Particles(sphere,[=](const TV& X){return TV(0.1,0);},0,
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>()+1.5);
            Add_Fixed_Corotated(1e2*scale_E,0.3);
        } break;
        case 3:{ // freefall circle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
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
                Seed_Particles(sphere,[=](const TV& X){return v0(s);},0,density,particles_per_cell);
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
            Seed_Particles(sphere,[=](const TV& X){return TV(0.5,0);},0,
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
            Seed_Particles(sphere1,[=](const TV& X){return TV(0.75,0);},0,density,particles_per_cell);
            SPHERE<TV> sphere2(TV(16,5),2);
            Seed_Particles(sphere2,[=](const TV& X){return TV(-0.75,0);},0,density,particles_per_cell);
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
                Seed_Particles(sphere,[=](const TV& X){return v0(s);},0,density,particles_per_cell);
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
            Seed_Particles(sphere,[=](const TV& X){return TV(3.0,0);},0,
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
            Seed_Particles(sphere,[=](const TV& X){return TV(3.0,0);},0,
                density,particles_per_cell);
            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            Add_Neo_Hookean(scale_E,0.425,&mpm_particles);

            SPHERE<TV> sphere2(TV(3.55,3.35),.3);
            TRIANGULATED_AREA<T>* ta=TESSELLATION::Generate_Triangles(sphere2,ceil(sphere2.radius/grid.dX.x*2));
            TRIANGULATED_AREA<T>& new_ta=Seed_Lagrangian_Particles(*ta,[=](const TV& X){return TV(-3.0,0);},
                0,density,true);
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
                Seed_Particles(sphere,[=](const TV& X){return TV(1,0);},0,
                    density,particles_per_cell);
                ARRAY<int> foo(IDENTITY_ARRAY<>(particles.number));
                Add_Fixed_Corotated(20*scale_E,0.3,&foo);}
            int N_sphere_particles=particles.number;
            {RANGE<TV> box(TV(.5,.2),TV(.6,.8));
                T density=1*scale_mass;
                Seed_Particles(box,0,0,density,particles_per_cell);
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
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            ARRAY<int> foo(IDENTITY_ARRAY<>(particles.number));
            Add_Fixed_Corotated(1*scale_E,0.3,&foo);
            Add_Gravity(TV(0,-9.8));
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1,false);
            MPM_COLLISION_IMPLICIT_OBJECT<TV>* bottom=dynamic_cast<MPM_COLLISION_IMPLICIT_OBJECT<TV>*>(collision_objects(3));
            bottom->func_frame=[this](T time){return FRAME<TV>(TV(0,(T).75-abs((T).15*time*scale_speed-(T).75)));};
            bottom->func_twist=[this](T time){return TWIST<TV>(TV(0,-sign((T).15*time*scale_speed-(T).75)*(T).15*scale_speed),typename TV::SPIN());};
        } break;
        case 12:{ // jello in a shrinking penalty box
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(TV(.3,.2),TV(.7,.4));
            T density=2*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            ARRAY<int> foo(IDENTITY_ARRAY<>(particles.number));
            Add_Fixed_Corotated(20*scale_E,0.4,&foo,false);
            Add_Gravity(TV(0,-9.8));
            Add_Walls(-1,COLLISION_TYPE::separate,1.9,.1,true);
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
            SEGMENTED_CURVE_2D<T>& new_sc=Seed_Lagrangian_Particles(*sc,[=](const TV& X){return TV(0.0,0);},0,density,true);
            SURFACE_TENSION_FORCE<TV>* stf=new SURFACE_TENSION_FORCE<TV>(new_sc,(T)0.01);
            Add_Force(*stf);
        } break;
        case 14:{ // surface tension static circle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            T density=2*scale_mass;
            SPHERE<TV> sphere(TV(.5,.5),.2);
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            bool no_mu=true;
            Add_Fixed_Corotated(scale_E*10,0.3,&mpm_particles,no_mu);
            use_surface_tension=true;
        } break;
        case 15:{ // surface tension become a circle
            //./mpm 15 -affine -resolution 64 -max_dt 1e-3 
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            T density=2*scale_mass;
            RANGE<TV> box(TV(.3,.11),TV(.5,.31));
            Seed_Particles(box,0,0,density,particles_per_cell);
            RANGE<TV> box2(TV(.6,.15),TV(.8,.25));
            Seed_Particles(box2,[=](const TV& X){return TV(-0.1,0);},0,density,particles_per_cell);
            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            bool no_mu=true;
            Add_Fixed_Corotated(scale_E*10,0.3,&mpm_particles,no_mu);
            use_surface_tension=true;
        } break;
        case 16:{ // oscillating circle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            use_oldroyd=true;
            this->inv_Wi=(T)0;
            particles.Store_S(use_oldroyd);            
            Seed_Particles(sphere,[=](const TV& X){return TV(0.1,0);},0,
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>()+1.5);
            particles.S.Fill(SYMMETRIC_MATRIX<T,2>()+sqr(1.5));
            OLDROYD_NEO_HOOKEAN<TV> *neo=new OLDROYD_NEO_HOOKEAN<TV>;
            neo->mu=38.462; // E=100, nu=0.3
            neo->lambda=57.692;
            MPM_OLDROYD_FINITE_ELEMENTS<TV> *fe=new MPM_OLDROYD_FINITE_ELEMENTS<TV>(force_helper,*neo,gather_scatter,0,this->inv_Wi,quad_F_coeff);
            Add_Force(*fe);
        } break;
        case 17:{ // spring test
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SEGMENTED_CURVE_2D<T>* sc=SEGMENTED_CURVE_2D<T>::Create();
            TV S(0.3,0.6);
            TV E(0.7,0.5);
            int N=20;
            TV D=(E-S)/(T)(N-1);
            T density=1*scale_mass;
            sc->particles.Add_Elements(N);
            sc->Update_Number_Nodes();
            for(int n=0;n<N-1;n++){
                TV X=S+D*n;
                sc->particles.X(n)=X;
                sc->mesh.elements.Append(TV_INT(n,n+1));}
            sc->particles.X(N-1)=E;
            SEGMENTED_CURVE_2D<T>& new_sc=Seed_Lagrangian_Particles(*sc,[=](const TV& X){return TV(0.0,0);},0,density,true);
            LINEAR_SPRINGS<TV>* stf=new LINEAR_SPRINGS<TV>(particles,new_sc.mesh);
            stf->Set_Restlength_From_Particles();
            stf->Set_Stiffness((T)10);
            stf->Set_Damping((T)0);
            Add_Force(*stf);
            for(int n=0;n<N;n++)
                for(int d=0;d<TV::m;d++)
                    particles.X(n)(d)+=random.Get_Uniform_Number(-0.2,0.2);
        } break;
        case 21:{ // circle with random initial velocities
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(-3,-3),TV(4,4)),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            random.Fill_Uniform(particles.V,-1,1);
            Add_Fixed_Corotated(1e3*scale_E,0.3);
        } break;
        case 22:{ // (fluid test) pool of water 
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(grid.dX*2,TV(1-2*grid.dX(0),0.25));
            T density=2*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Fixed_Corotated(1e3*scale_E,0.3);
            Add_Gravity(TV(0,-1.8));
            if(this->kkt) particles.lambda.Fill(FLT_MAX);
            Add_Fluid_Wall(new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(RANGE<TV>(TV(-5,-5),TV(5,0.1))));
        } break;
        case 23:{ // (fluid test) dam break 
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(grid.dX*2,TV(0.2,0.75));
            T density=2*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Fixed_Corotated(1e3*scale_E,0.3);
            Add_Gravity(TV(0,-1.8));
        } break;
        case 24:{ // (fluid test) circle drop 
            // one: ./mpm -kkt -scale_E 0
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.7),.2);
            T density=2*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            Add_Fixed_Corotated(1e3*scale_E,0.3);
            Add_Gravity(TV(0,-1.8));
            if(this->kkt) particles.lambda.Fill(FLT_MAX);
        } break;
        case 25:{ // (fluid test) pool of water w/ single particle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(grid.dX*(T).5,TV(1-grid.dX(0)*(T).5,0.25));
            T density=2*scale_mass;
            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Particle(TV(.5,.9),0,0,mass,volume);
            Add_Gravity(TV(0,-1.8));
        } break;
        case 26:{ // Rayleigh Taylor
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(grid.dX*2,TV(1-2*grid.dX(0),0.20));
            T density=2*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            density*=10;
            box+=TV(0,0.20-2*grid.dX(0));
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Gravity(TV(0,-1.8));
        } break;
        case 98:{ // full box
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(grid.dX*(T)2,TV::All_Ones_Vector()-grid.dX*2);
            T density=2*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Gravity(TV(0,-1.8));
        } break;
        case 99:{ // single particle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            T density=2*scale_mass;
            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            Add_Particle(TV(.8,.5),0,0,mass,volume);
            Add_Gravity(TV(0,-1.8));
        } break;
        case 27:{ // drop an oldroyd-b to a ground
            grid.Initialize(TV_INT(resolution*2,resolution),RANGE<TV>(TV(-1,0),TV(1,1)),true);
            Add_Collision_Object(RANGE<TV>(TV(-5,-5),TV(5,.1)),COLLISION_TYPE::stick,0);
            SPHERE<TV> sphere(TV(.5,.5),.2);
            T density=2*scale_mass;
            use_oldroyd=true;
            this->inv_Wi=(T)100;
            particles.Store_S(use_oldroyd);            
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>()+1);particles.S.Fill(SYMMETRIC_MATRIX<T,2>()+sqr(1));
            LOG::cout<<particles.F<<std::endl<<std::endl;;
            LOG::cout<<particles.S<<std::endl;
            OLDROYD_NEO_HOOKEAN<TV> *neo=new OLDROYD_NEO_HOOKEAN<TV>;
            neo->mu=38.462; // E=100, nu=0.3
            neo->lambda=57.692;
            MPM_OLDROYD_FINITE_ELEMENTS<TV> *fe=new MPM_OLDROYD_FINITE_ELEMENTS<TV>(force_helper,*neo,gather_scatter,0,this->inv_Wi,quad_F_coeff);
            Add_Force(*fe);
            Add_Gravity(TV(0,-1.8));
        } break;
        case 29:{ // drop an oldroyd-b to a ground SCA energy
            grid.Initialize(TV_INT(resolution*2,resolution),RANGE<TV>(TV(-1,0),TV(1,1)),true);
            Add_Collision_Object(RANGE<TV>(TV(-5,-5),TV(5,.1)),COLLISION_TYPE::stick,0);
            SPHERE<TV> sphere(TV(.5,.5),.2);
            T density=2*scale_mass;
            use_oldroyd=true;
            this->inv_Wi=(T)0;
            particles.Store_S(use_oldroyd);            
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>()+1);particles.S.Fill(SYMMETRIC_MATRIX<T,2>()+sqr(1));
            VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV> *neo=new VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV>;
            neo->mu=0;
            neo->lambda=57.692;
            MPM_OLDROYD_FINITE_ELEMENTS<TV> *fe=new MPM_OLDROYD_FINITE_ELEMENTS<TV>(force_helper,*neo,gather_scatter,0,this->inv_Wi,quad_F_coeff);
            Add_Force(*fe);
            Add_Force(*new MPM_VISCOSITY<TV>(force_helper,gather_scatter,0,0.001));
            Add_Gravity(TV(0,-1.8));
        } break;
        case 28:{ // newton convergence problem: ./mpm 28 -affine -max_dt 1e-3 | grep converge
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            Add_Collision_Object(RANGE<TV>(TV(-5,-5),TV(5,.1)),COLLISION_TYPE::slip,0);
            Add_Collision_Object(RANGE<TV>(TV(-5,-5),TV(.1,5)),COLLISION_TYPE::slip,0);
            RANGE<TV> box(TV(.21,.21),TV(.6,.5));
            T density=2*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Gravity(TV(0,-1.8));
            Add_Fixed_Corotated(1.71*scale_E,0.4);
        } break;
        case 30:{ // pinned rotating circle
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            VECTOR<T,1> angular_velocity(0.4*scale_speed);
            T density=2*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,penalty_collisions_stiffness,
                penalty_damping_stiffness);
            T pin_radius=sphere.radius*(T).8;
            for(int i=0;i<particles.X.m;i++)
                if((particles.X(i)-sphere.center).Magnitude_Squared()>sqr(pin_radius)){
                    TV dx=particles.X(i)-sphere.center;
                    pinning_force->Add_Target(i,
                        [=](T time){
                            ROTATION<TV> rot=ROTATION<TV>::From_Rotation_Vector(angular_velocity*time);
                            return rot.Rotate(dx)+sphere.center;});}
            Add_Force(*pinning_force);
            Add_Neo_Hookean(1e3*scale_E,0.3);
        } break;
        case 32:{ // colliding
            grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(-1,-1),TV(2,2)),true);
            T d_small_1=1*scale_mass,d_small_2=1*scale_mass,d_large_1=d_small_1*.9,d_large_2=d_small_2*.9;
            SPHERE<TV> large1(TV(.2,.5),.1);
            Seed_Particles(large1,[=](const TV& X){return TV(0.75,0);},0,d_large_1,particles_per_cell);
            SPHERE<TV> large2(TV(.6,.5),.1);
            Seed_Particles(large2,[=](const TV& X){return TV();},0,d_large_2,particles_per_cell);
            for(int k=0;k<particles.number;k++)
                if((particles.X(k)-large1.center).Magnitude_Squared()<sqr(large1.radius*0.6)
                    || (particles.X(k)-large2.center).Magnitude_Squared()<sqr(large2.radius*0.6))
                    particles.deletion_list.Append(k);
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
            for(int i=0;i<particles.X.m;i++) (*color_attribute)(i)=VECTOR<T,3>(0,1,1);
            particles.Delete_Elements_On_Deletion_List();
            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            Add_Neo_Hookean(2*scale_E,0.425,&mpm_particles);
            int old_pn=particles.number;
            SPHERE<TV> small1(large1.center,large1.radius*.6);
            Seed_Particles(small1,[=](const TV& X){return TV(0.75,0);},0,d_small_1,particles_per_cell);
            SPHERE<TV> small2(large2.center,large2.radius*.6);
            Seed_Particles(small2,[=](const TV& X){return TV();},0,d_small_2,particles_per_cell);
            ARRAY<int> mpm_particles2;for(int i=old_pn;i<particles.number;i++) mpm_particles2.Append(i);
            Add_Neo_Hookean(2000*scale_E,0.3,&mpm_particles2);
            for(int i=old_pn;i<particles.X.m;i++) (*color_attribute)(i)=VECTOR<T,3>(1,1,1);
        } break;
        case 33:{ // sand box drop
            particles.Store_Fp(true);
            particles.Store_Lame(true);

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1)),COLLISION_TYPE::separate,10);

            T density=(T)1281*scale_mass;
            T E=5000*scale_E,nu=.4;
            if(!use_theta_c) theta_c=0.01;
            if(!use_theta_s) theta_s=.00001;
            if(!use_hardening_factor) hardening_factor=80;
            if(!use_max_hardening) max_hardening=5;
            Add_Clamped_Plasticity(*new COROTATED_FIXED<T,TV::m>(E,nu),theta_c,theta_s,max_hardening,hardening_factor,0);
            RANGE<TV> box(TV(.45,.11),TV(.55,.31));
            Seed_Particles(box,0,0,density,particles_per_cell);
            T mu=E/(2*(1+nu));
            T lambda=E*nu/((1+nu)*(1-2*nu));
            particles.mu.Fill(mu);
            particles.mu0.Fill(mu);
            particles.lambda.Fill(lambda);
            particles.lambda0.Fill(lambda);
            Add_Gravity(TV(0,-9.8));
        } break;
        case 34:{ // drip drop
            grid.Initialize(TV_INT(1,2)*resolution,RANGE<TV>(TV(0,-1),TV(1,1)),true);
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1,true);
            Add_Gravity(TV(0,-1));
            T density=2*scale_mass;
            RANGE<TV> box(TV(.4,.5),TV(.6,.85));
            Seed_Particles(box,[=](const TV& X){return TV(-0.0,0);},0,density,particles_per_cell);
            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            bool no_mu=true;
            Add_Fixed_Corotated(scale_E*20,0.3,&mpm_particles,no_mu);
            Add_Force(*new MPM_VISCOSITY<TV>(force_helper,gather_scatter,0,0.000001));
            use_surface_tension=true;
            PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,penalty_collisions_stiffness,penalty_damping_stiffness);
            for(int i=0;i<particles.X.m;i++)
                if(particles.X(i)(1)>0.8){TV x=particles.X(i);pinning_force->Add_Target(i,[=](T time){return x;});}
            Add_Force(*pinning_force);
        } break;
        case 35:{ // snow wedge
            // ./mpm 35 -flip 0.95 -max_dt .005 -cfl .1 -resolution 200
            particles.Store_Fp(true);
            particles.Store_Lame(true);

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            ORIENTED_BOX<TV> wedge(RANGE<TV>(TV(),TV(0.2,0.2)),ROTATION<TV>::From_Angle(0.25*M_PI),TV(0.5,0.4-sqrt(2.0)*0.1));
            RANGE<TV> ground(TV(-1,0),TV(2,0.1));
            if(use_penalty_collisions){
                Add_Penalty_Collision_Object(wedge);
                Add_Penalty_Collision_Object(ground);}
            else{
                Add_Collision_Object(ground,COLLISION_TYPE::separate,1);
                Add_Collision_Object(new ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<TV> >(wedge),COLLISION_TYPE::separate,1);}

            T density=(T)2*scale_mass;
            int number_of_particles=20000;
            T E=40*scale_E,nu=.2;
            if(!use_theta_c) theta_c=0.015;
            if(!use_theta_s) theta_s=.005;
            if(!use_hardening_factor) hardening_factor=7;
            if(!use_max_hardening) max_hardening=FLT_MAX;
            Add_Clamped_Plasticity(*new COROTATED_FIXED<T,TV::m>(E,nu),theta_c,theta_s,max_hardening,hardening_factor,0);
            RANGE<TV> box(TV(.3,.7),TV(.7,.9));
            std::ifstream ifs("particles.dat");
            if(ifs.is_open()){
                std::string foobar;
                ifs>>foobar>>foobar>>foobar>>foobar;
                T x,y,mass,volume;
                while(ifs>>x>>y>>mass>>volume && !ifs.eof())
                    Add_Particle(TV(x,y),0,0,mass,volume);
                PHYSBAM_ASSERT(particles.number==number_of_particles);}
            else{
                PHYSBAM_WARNING("Couldn't open 'particles.dat'. Falling back to using random particle positions.");
                Seed_Particles(box,0,0,density,number_of_particles*grid.dX.Product()/box.Size());}
            T mu=E/(2*(1+nu));
            T lambda=E*nu/((1+nu)*(1-2*nu));
            particles.mu0.Fill(mu);
            particles.mu.Fill(mu);
            particles.lambda0.Fill(lambda);
            particles.lambda.Fill(lambda);
            Add_Gravity(TV(0,-2));
        } break;
        case 36:{ // split
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV(0,0);},0,
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>(1,0,0,1));
            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            Add_Fixed_Corotated(scale_E*10,0.3,&mpm_particles);
            PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,penalty_collisions_stiffness,penalty_damping_stiffness);
            for(int i=0;i<particles.X.m;i++)
                if(particles.X(i)(0)<0.25 || particles.X(i)(0)>90.75){TV x=particles.X(i);pinning_force->Add_Target(i,[=](T time){return x;});}
            Add_Force(*pinning_force);
            SEGMENTED_CURVE_2D<T>* sc=SEGMENTED_CURVE_2D<T>::Create();
            sc->particles.Add_Elements(particles.number);
            sc->Update_Number_Nodes();
            for(int p=0;p<particles.number;p++){
                sc->particles.X(p)=particles.X(p);
                for(int q=0;q<particles.number;q++){
                    TV L=particles.X(p),R=particles.X(q);
                    if(L(0)<=0.5 && R(0)>0.5&&  L(1)>0.5 && R(1)>0.5 && (L-R).Magnitude_Squared()<sqr(grid.dX.Min()*2)){
                        sc->mesh.elements.Append(TV_INT(p,q));}}}
            LINEAR_SPRINGS<TV>* stf=new LINEAR_SPRINGS<TV>(particles,sc->mesh);
            ARRAY<T> r(sc->mesh.elements.m);r.Fill(grid.dX.Min()*2.1);
            stf->Set_Restlength(r);
            stf->Set_Stiffness((T)1);
            stf->Set_Damping((T)0.1);
            Add_Force(*stf);
        } break;
        case 37:{ // sand box drop, better paramaters, with Hencky, usage: mpm 38 -resolution 100 -plastic_newton_iterations 100 -plastic_newton_tolerance 1e-8
            particles.Store_Fp(true);
            particles.Store_Lame(true);

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            if(use_penalty_collisions)
                Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1)),0.9);
            else
                Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1)),COLLISION_TYPE::stick,0);

            T density=(T)1281*scale_mass;
            T E=35.37e6*scale_E,nu=.4;
            Add_St_Venant_Kirchhoff_Hencky_Strain(E,nu);
//            this->plasticity=new MPM_DRUCKER_PRAGER<TV>(friction_angle,cohesion);
            RANGE<TV> box(TV(.4,.1001),TV(.45,.6001));
            Seed_Particles(box,0,0,density,particles_per_cell);
            T mu=E/(2*(1+nu));
            T lambda=E*nu/((1+nu)*(1-2*nu));
            particles.mu.Fill(mu);
            particles.mu0.Fill(mu);
            particles.lambda.Fill(lambda);
            particles.lambda0.Fill(lambda);
            Add_Gravity(TV(0,-9.81));
        } break;
        case 38:{ // sand box drop, wide
            particles.Store_Fp(true);
            particles.Store_Lame(true);

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            if(use_penalty_collisions)
                Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1)),0.9);
            else
                Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1)),COLLISION_TYPE::stick,0);

            T density=(T)1281*scale_mass;
            T E=35.37e6*scale_E,nu=.4;
            Add_St_Venant_Kirchhoff_Hencky_Strain(E,nu);
//            this->plasticity=new MPM_DRUCKER_PRAGER<TV>(friction_angle,cohesion);
            T l0=0.05;
            T h0=l0*8;
            T gap=grid.dX(1)*0.01;
            RANGE<TV> box(TV(.5-l0,.1+gap),TV(.5+l0,.1+gap+h0));
            Seed_Particles(box,0,0,density,particles_per_cell);
            T mu=E/(2*(1+nu));
            T lambda=E*nu/((1+nu)*(1-2*nu));
            particles.mu.Fill(mu);
            particles.mu0.Fill(mu);
            particles.lambda.Fill(lambda);
            particles.lambda0.Fill(lambda);
            Add_Gravity(TV(0,-9.81));
        } break;
        case 39:{ // DP on wedge
            particles.Store_Fp(true);
            particles.Store_Lame(true);

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            ORIENTED_BOX<TV> wedge(RANGE<TV>(TV(),TV(0.2,0.2)),ROTATION<TV>::From_Angle(0.25*M_PI),TV(0.5,0.4-sqrt(2.0)*0.1));
            RANGE<TV> ground(TV(-1,0),TV(2,0.1));
            if(use_penalty_collisions){
                Add_Penalty_Collision_Object(wedge);
                Add_Penalty_Collision_Object(ground);}
            else{
                Add_Collision_Object(ground,COLLISION_TYPE::separate,1);
                Add_Collision_Object(new ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<TV> >(wedge),COLLISION_TYPE::separate,1);}

            T density=(T)2*scale_mass;
            int number_of_particles=20000;
            T E=40*scale_E,nu=.2;
            Add_St_Venant_Kirchhoff_Hencky_Strain(E,nu);
//            this->plasticity=new MPM_DRUCKER_PRAGER<TV>(friction_angle,cohesion);
            RANGE<TV> box(TV(.3,.7),TV(.7,.9));
            std::ifstream ifs("particles.dat");
            if(ifs.is_open()){
                std::string foobar;
                ifs>>foobar>>foobar>>foobar>>foobar;
                T x,y,mass,volume;
                while(ifs>>x>>y>>mass>>volume && !ifs.eof())
                    Add_Particle(TV(x,y),0,0,mass,volume);
                PHYSBAM_ASSERT(particles.number==number_of_particles);}
            else{
                PHYSBAM_WARNING("Couldn't open 'particles.dat'. Falling back to using random particle positions.");
                Seed_Particles(box,0,0,density,number_of_particles*grid.dX.Product()/box.Size());}
            T mu=E/(2*(1+nu));
            T lambda=E*nu/((1+nu)*(1-2*nu));
            particles.mu0.Fill(mu);
            particles.mu.Fill(mu);
            particles.lambda0.Fill(lambda);
            particles.lambda.Fill(lambda);
            Add_Gravity(TV(0,-2));
        } break;
        case 40:
        case 41:
        case 42:
        case 43:
        case 44:
        case 45:
        case 46:
        case 47:
        case 48:
        case 49:{ // Mast paper
            particles.Store_Fp(true);
            particles.Store_Lame(true);

            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            if(use_penalty_collisions)
                Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1)),0.9);
            else
                Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1)),COLLISION_TYPE::stick,0);

            T density=(T)2200*scale_mass;
            T E=35.37e6*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Add_Drucker_Prager_Case(E,nu,test_number-40);
            T l0=0.05;
            T h0=l0*8;
            T gap=grid.dX(1)*0.01;
            RANGE<TV> box(TV(.5-l0,.1+gap),TV(.5+l0,.1+gap+h0));
            Seed_Particles(box,0,0,density,particles_per_cell);
            T mu=E/(2*(1+nu));
            T lambda=E*nu/((1+nu)*(1-2*nu));
            particles.mu.Fill(mu);
            particles.mu0.Fill(mu);
            particles.lambda.Fill(lambda);
            particles.lambda0.Fill(lambda);
            Add_Gravity(TV(0,-9.81));
        } break;
        case 50:
        case 51:{ //lambda particles
            //usage:./mpm 50 -use_exp_F -max_dt 1e-3 -resolution 100
            particles.Store_Fp(true);
            particles.Store_Lame(true);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            //particles.Store_Plastic_Deformation(true);
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box(),true);
            if(use_penalty_collisions){
                Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1)),0.9);
                Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(.1,1.5)),0.9);
                Add_Penalty_Collision_Object(RANGE<TV>(TV(.9,-1),TV(1.5,1.5)),0.9);}
            else{
                Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1)),COLLISION_TYPE::stick,0);
                Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(.1,1.5)),COLLISION_TYPE::stick,0);
                Add_Collision_Object(RANGE<TV>(TV(.9,-1),TV(1.5,1.5)),COLLISION_TYPE::stick,0);}
            
            T density=(T)2200*scale_mass;
            T E=35.37e6*scale_E,nu=.3;
            T mu=E/(2*(1+nu));
            T lambda=E*nu/((1+nu)*(1-2*nu));
            //this->plasticity=new MPM_DRUCKER_PRAGER_HARDENING<TV>(35,0,0,0);
            T gap=grid.dX(1)*0.01;
            RANGE<TV> box(TV(.1+gap,.1+gap),TV(.3,.75));
            Seed_Particles(box,0,0,density,particles_per_cell);
            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager(E,nu,(T)35,&sand_particles);
            particles.mu.Fill(mu);
            particles.mu0.Fill(mu);
            particles.lambda.Fill(lambda);
            particles.lambda0.Fill(lambda);
            
            if(test_number==51){
                T El=5000*scale_E,nul=.3;
                Add_Lambda_Particles(&sand_particles,El,nul,true);}

            Add_Gravity(TV(0,-9.81));
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
    if(test_number==12){
        if(time>=10/24.0){
            lagrangian_forces.Delete_Pointers_And_Clean_Memory();
            this->deformable_body_collection.structures.Delete_Pointers_And_Clean_Memory();
            this->output_structures_each_frame=true;
            Add_Walls(-1,COLLISION_TYPE::separate,1.9,.1+(T)(time-10/24.0)*0.08,true);
            Add_Gravity(TV(0,-9.8));}}
    if(test_number==36){
        delete lagrangian_forces(lagrangian_forces.m-1);
        lagrangian_forces.Remove_End();
        SEGMENTED_CURVE_2D<T>* sc=SEGMENTED_CURVE_2D<T>::Create();
        sc->particles.Add_Elements(particles.number);
        sc->Update_Number_Nodes();
        for(int p=0;p<particles.number;p++){
            sc->particles.X(p)=particles.X(p);
            for(int q=0;q<particles.number;q++){
                TV L=particles.X(p),R=particles.X(q);
                if(L(0)<=0.5 && R(0)>0.5 && L(1)>0.5 && R(1)>0.5 && (L-R).Magnitude_Squared()<sqr(grid.dX.Min()*2)){
                    sc->mesh.elements.Append(TV_INT(p,q));}}}
        LINEAR_SPRINGS<TV>* stf=new LINEAR_SPRINGS<TV>(particles,sc->mesh);
        ARRAY<T> r(sc->mesh.elements.m);r.Fill(grid.dX.Min()*2.1);
        stf->Set_Restlength(r);
        stf->Set_Stiffness((T)1);
        stf->Set_Damping((T)0);
        Add_Force(*stf);}
    if(use_surface_tension){

        bool use_bruteforce=true;
        bool use_kdtree=false;

        // Remove old surface particles
        int N_non_surface=particles.number-Nsurface;
        for(int k=N_non_surface;k<particles.number;k++){
            particles.Add_To_Deletion_List(k);
            int m=steal(k-N_non_surface);
            TV old_momentum=particles.mass(m)*particles.V(m)+particles.mass(k)*particles.V(k);
            particles.mass(m)+=particles.mass(k);
            particles.volume(m)+=particles.volume(k);
            particles.V(m)=old_momentum/particles.mass(m);
            if(use_affine) particles.B(m)=(particles.B(m)+particles.B(k))*0.5;}
        LOG::cout<<"deleting "<<Nsurface<<" particles..."<<std::endl;
        particles.Delete_Elements_On_Deletion_List();
        lagrangian_forces.Delete_Pointers_And_Clean_Memory();
        this->deformable_body_collection.structures.Delete_Pointers_And_Clean_Memory();

        // Dirty hack of rasterizing mass
        this->simulated_particles.Remove_All();
        for(int p=0;p<this->particles.number;p++)
            if(this->particles.valid(p))
                this->simulated_particles.Append(p);
        this->particle_is_simulated.Remove_All();
        this->particle_is_simulated.Resize(this->particles.X.m);
        this->particle_is_simulated.Subset(this->simulated_particles).Fill(true);
        this->weights->Update(this->particles.X);
        this->gather_scatter.Prepare_Scatter(this->particles);
        MPM_PARTICLES<TV>& my_particles=this->particles;
#pragma omp parallel for
        for(int i=0;i<this->mass.array.m;i++)
            this->mass.array(i)=0;
        this->gather_scatter.template Scatter<int>(false,0,
            [this,&my_particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
            {
                T w=it.Weight();
                TV_INT index=it.Index();
                this->mass(index)+=w*my_particles.mass(p);
            });

        // Marching cube
        SEGMENTED_CURVE_2D<T>* surface=SEGMENTED_CURVE_2D<T>::Create();
        MARCHING_CUBES<TV>::Create_Surface(*surface,grid,mass,particles.mass(0)*.7);

        // Improve surface quality
        T min_edge_length=FLT_MAX;
        for(int k=0;k<surface->mesh.elements.m;k++){
            int a=surface->mesh.elements(k)(0),b=surface->mesh.elements(k)(1);
            TV A=surface->particles.X(a),B=surface->particles.X(b);
            T l2=(A-B).Magnitude_Squared();
            if(l2<min_edge_length) min_edge_length=l2;}
        min_edge_length=sqrt(min_edge_length);
        LOG::cout<<"Marching cube min edge length: "<<min_edge_length<<std::endl;

        // Seed surface particles
        int Nold=particles.number;
        Nsurface=surface->particles.number;
        LOG::cout<<"adding "<<Nsurface<<" particles..."<<std::endl;
        SEGMENTED_CURVE_2D<T>& new_sc=Seed_Lagrangian_Particles(*surface,0,0,(T)0.001,true);
        ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
        for(int i=Nold;i<particles.X.m;i++) (*color_attribute)(i)=VECTOR<T,3>(0,1,0);

        // Build K-d tree for non-surface particles
        LOG::cout<<"building kdtree..."<<std::endl;
        KD_TREE<TV> kdtree;
        ARRAY<TV> nodes(Nold);
        if(use_kdtree){
            for(int p=0;p<Nold;p++) nodes(p)=particles.X(p);
            kdtree.Create_Left_Balanced_KD_Tree(nodes);}

        // Assign physical quantities
        steal.Clean_Memory();
        for(int k=Nold;k<particles.number;k++){
            T dist2=FLT_MAX;
            int m=-1;

            // Find closest interior particle using brute force
            if(use_bruteforce){
                for(int q=0;q<Nold;q++){
                    T dd=(particles.X(q)-particles.X(k)).Magnitude_Squared();
                    if(dd<dist2){dist2=dd;m=q;}}}

            // Find closest interior particle using kdtree
            int number_of_points_in_estimate=1;
            ARRAY<int> points_found(number_of_points_in_estimate);
            ARRAY<T> squared_distance_to_points_found(number_of_points_in_estimate);
            if(use_kdtree){
                int number_of_points_found;T max_squared_distance_to_points_found;
                kdtree.Locate_Nearest_Neighbors(particles.X(k),FLT_MAX,points_found,
                    squared_distance_to_points_found,number_of_points_found,max_squared_distance_to_points_found,nodes);
                PHYSBAM_ASSERT(number_of_points_found==number_of_points_in_estimate);}

            // Debug kdtree
            if(use_bruteforce && use_kdtree && m!=points_found(0)){
                LOG::cout<<"Disagree!"<<std::endl; PHYSBAM_FATAL_ERROR();}

            if(use_kdtree) m=points_found(0);

            T split_mass=particles.mass(m)*.5;
            T split_volume=particles.volume(m)*.5;
            TV com=particles.X(m);
            particles.mass(m)=split_mass;
            particles.volume(m)=split_volume;
            particles.mass(k)=split_mass;
            particles.volume(k)=split_volume;
            particles.V(k)=particles.V(m);
            if(use_affine) particles.B(k)=particles.B(m);
            steal.Append(m);}

        // Add surface tension force
        SURFACE_TENSION_FORCE<TV>* stf=new SURFACE_TENSION_FORCE<TV>(new_sc,(T)1e-2);
        Add_Force(*stf);
    }
    
    if(test_number==34){
        Add_Walls(-1,COLLISION_TYPE::separate,.3,.1,true);
        Add_Gravity(TV(0,-1));
        PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,penalty_collisions_stiffness,penalty_damping_stiffness);
        for(int i=0;i<particles.X.m;i++)
            if(particles.X(i)(1)>0.8){TV x=particles.X(i);pinning_force->Add_Target(i,[=](T time){return x;});}
        Add_Force(*pinning_force);
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
