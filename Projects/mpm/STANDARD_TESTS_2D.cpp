//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/KD_TREE.h>
#include <Core/Math_Tools/pow.h>
#include <Core/Math_Tools/RANGE_ITERATOR.h>
#include <Core/Matrices/FRAME.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/HOURGLASS.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INTERSECTION.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INVERT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_UNION.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_UTILITIES.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Seeding/POISSON_DISK.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/PINNING_FORCE.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Constitutive_Models/QUASI_INCOMPRESSIBLE_FORCE.h>
#include <Deformables/Constitutive_Models/TAIT_PRESSURE_FORCE.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_DRUCKER_PRAGER.h>
#include <Hybrid_Methods/Forces/MPM_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_GRAVITY.h>
#include <Hybrid_Methods/Forces/MPM_OLDROYD_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_VISCOSITY.h>
#include <Hybrid_Methods/Forces/OLDROYD_NEO_HOOKEAN.h>
#include <Hybrid_Methods/Forces/VOLUME_PRESERVING_OB_NEO_HOOKEAN.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include "STANDARD_TESTS_2D.h"
#include <fstream>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS<VECTOR<T,2> >::
STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_BASE<TV>(stream_type_input,parse_args),
    foo_T1(0),foo_T2(0),foo_T3(0),foo_T4(0),
    use_foo_T1(false),use_foo_T2(false),use_foo_T3(false),use_foo_T4(false)
{
    parse_args.Add("-fooT1",&foo_T1,&use_foo_T1,"T1","a scalar");
    parse_args.Add("-fooT2",&foo_T2,&use_foo_T2,"T2","a scalar");
    parse_args.Add("-fooT3",&foo_T3,&use_foo_T3,"T3","a scalar");
    parse_args.Add("-fooT4",&foo_T4,&use_foo_T4,"T4","a scalar");
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
    switch(test_number)
    {
        case 1:{ // rotating circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            VECTOR<T,1> angular_velocity(0.4*scale_speed/s);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return angular_velocity.Cross(X-sphere.center);},
                [=](const TV&){return MATRIX<T,2>::Cross_Product_Matrix(angular_velocity);},
                density,particles_per_cell);
            T total_mass=particles.mass.Sum();
            TV total_momentum=particles.V.Weighted_Sum(particles.mass);
            TV dV=total_momentum/total_mass;
            particles.V-=dV;
            T E=1e3*unit_p*scale_E;
            T nu=0.3;
            Add_Neo_Hookean(E,nu);
        } break;
        case 2:{ // oscillating circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV(0.1,0)*(m/s);},0,
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>()+1.5);
            Add_Fixed_Corotated(1e2*unit_p*scale_E,0.3);
        } break;
        case 3:{ // freefall circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-9.8));
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1*m,false);
//            Add_Fixed_Corotated(1e2*unit_p*scale_E,0.3);
        } break;
        case 4:{ // colliding of two rings
            Set_Grid(RANGE<TV>(TV(),TV(0.48,0.48))*m,TV_INT()+1,48);
            ARRAY<SPHERE<TV> > spheres;
            ARRAY<TV> v0;
            ARRAY<T> r;
            spheres.Append(SPHERE<TV>(TV(0.1,0.24)*m,0.04*m));
            spheres.Append(SPHERE<TV>(TV(0.4,0.24)*m,0.04*m));
            v0.Append(TV(50,0)*(m/s));
            v0.Append(TV(-50,0)*(m/s));
            r.Append(0.03*m);
            r.Append(0.03*m);
            for(int s=0;s<spheres.m;s++){
                SPHERE<TV>& sphere=spheres(s);
                T density=1010*unit_rho*scale_mass;
                int last=particles.number;
                Seed_Particles(sphere,[=](const TV& X){return v0(s);},0,density,particles_per_cell);
                for(int k=last;k<particles.number;k++){
                    if((particles.X(k)-sphere.center).Magnitude_Squared()<sqr(r(s))){
                        particles.deletion_list.Append(k);}}}
            particles.Delete_Elements_On_Deletion_List();
            Add_Neo_Hookean(0.073e9*unit_p*scale_E,0.4);
        } break;
        case 5:{ // rebound of an elastic cylinder
            if(!user_resolution) resolution=10;
            T dx=(T)5/resolution*m;
            Set_Grid(RANGE<TV>(TV(),TV(15,5)*m).Thickened(dx*(T)4.5),TV_INT(3,1),TV_INT()+9,1,10);
            Add_Collision_Object(RANGE<TV>(TV(-5,-5)*m,TV(0+dx/2,15*m)),COLLISION_TYPE::slip,0);
            Add_Collision_Object(RANGE<TV>(TV(15*m-dx/2,-5*m),TV(20,15)*m),COLLISION_TYPE::slip,0);
            SPHERE<TV> sphere(TV(2.5,2.5)*m,1.5*m);
            T density=4*unit_rho*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV(0.5*m/s,0);},0,
                density,particles_per_cell);
            Add_Neo_Hookean(85.5*unit_p*scale_E,0.425); //solve({E/(2*(1+r))=30,E*r/((1+r)*(1-2*r))=170},{E,r});
        } break;
        case 6:{ // skew impact of two elastic cylinders
            if(!user_resolution) resolution=12;
            int third_resolution=(resolution+2)/3;
            T dx=(T)4/third_resolution;
            Set_Grid(RANGE<TV>(TV(),TV(20,12)*m).Thickened(dx*(T)4.5),TV_INT(5,3),TV_INT()+9,3,12);
            T density=5*unit_rho*scale_mass;
            SPHERE<TV> sphere1(TV(3,3)*m,2*m);
            Seed_Particles(sphere1,[=](const TV& X){return TV(0.75*(m/s),0);},0,density,particles_per_cell);
            SPHERE<TV> sphere2(TV(16,5)*m,2*m);
            Seed_Particles(sphere2,[=](const TV& X){return TV(-0.75*(m/s),0);},0,density,particles_per_cell);
            Add_Neo_Hookean(31.685*unit_p*scale_E,0.44022); //solve({E/(2*(1+r))=11,E*r/((1+r)*(1-2*r))=81},{E,r});
        } break;
        case 7:{ // ping-pong ring
            // ./mpm 7 -flip 0  -affine -midpoint -max_dt 1e-3 -cfl .1 -framerate 2400 -newton_tolerance 1e-5 -solver_tolerance 1e-5  -last_frame 240 -order 2 -print_stats | grep 'total'
            Set_Grid(RANGE<TV>(TV(),TV(0.48,0.48)*m),TV_INT()+1,480);
            Add_Collision_Object(RANGE<TV>(TV(-5,-5),TV(0.11,15))*m,COLLISION_TYPE::separate,0);
            Add_Collision_Object(RANGE<TV>(TV(0.3,-5),TV(20,15))*m,COLLISION_TYPE::separate,0);
            ARRAY<SPHERE<TV> > spheres;
            ARRAY<TV> v0;
            ARRAY<T> r;
            spheres.Append(SPHERE<TV>(TV(0.2,0.24)*m,0.04*m));
            v0.Append(TV(50*(m/s),0));
            r.Append(0.03*m);
            for(int s=0;s<spheres.m;s++){
                SPHERE<TV>& sphere=spheres(s);
                T density=1010*unit_rho*scale_mass;
                int last=particles.number;
                Seed_Particles(sphere,[=](const TV& X){return v0(s);},0,density,particles_per_cell);
                for(int k=last;k<particles.number;k++){
                    if((particles.X(k)-sphere.center).Magnitude_Squared()<sqr(r(s))){
                        particles.deletion_list.Append(k);}}}
            particles.Delete_Elements_On_Deletion_List();
            Add_Neo_Hookean(0.073e9*unit_p*scale_E,0.4);
        } break;
        case 8:{ // collision an elastic cylinder (TODO: fix description.)
            Set_Grid(RANGE<TV>(TV(),TV(5,5))*m,TV_INT()+1,TV_INT()+9,1,32);
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1*m,false);
            Add_Collision_Object(SPHERE<TV>(TV(4,3)*m,1*m),COLLISION_TYPE::separate,.3);
            SPHERE<TV> sphere(TV(2.55,3.55)*m,.3*m);
            T density=4*unit_rho*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV(3.0*(m/s),0);},0,
                density,particles_per_cell);
            Add_Neo_Hookean(unit_p*scale_E,0.425);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 9:{ // collision an elastic cylinder (TODO: fix description.)
            Set_Grid(RANGE<TV>(TV(),TV(5,5))*m,TV_INT()+1,TV_INT()+9,1,32);
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1*m,false);
            SPHERE<TV> sphere(TV(2.55,3.55)*m,.3*m);
            T density=4*unit_rho*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV(3.0*(m/s),0);},0,
                density,particles_per_cell);
            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            Add_Neo_Hookean(unit_p*scale_E,0.425,&mpm_particles);

            SPHERE<TV> sphere2(TV(3.55,3.35)*m,.3*m);
            TRIANGULATED_AREA<T>* ta=TESSELLATION::Generate_Triangles(sphere2,ceil(sphere2.radius/grid.dX.x*2));
            TRIANGULATED_AREA<T>& new_ta=Seed_Lagrangian_Particles(*ta,[=](const TV& X){return TV(-3.0*(m/s),0);},
                0,density,true);
            Add_Neo_Hookean(new_ta,unit_p*scale_E,0.425);

            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 10:{ // mpm projectile vs end-holded wall
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1*m,false);
            Add_Collision_Object(RANGE<TV>(TV(.45,.75),TV(.65,.85))*m,COLLISION_TYPE::stick,0);
            Add_Collision_Object(RANGE<TV>(TV(.45,.15),TV(.65,.25))*m,COLLISION_TYPE::stick,0);
            {SPHERE<TV> sphere(TV(.2,.5)*m,.06*m);
                T density=2*unit_rho*scale_mass;
                Seed_Particles(sphere,[=](const TV& X){return TV(1,0)*(m/s);},0,
                    density,particles_per_cell);
                ARRAY<int> foo(IDENTITY_ARRAY<>(particles.number));
                Add_Fixed_Corotated(20*unit_p*scale_E,0.3,&foo);}
            int N_sphere_particles=particles.number;
            {RANGE<TV> box(TV(.5,.2)*m,TV(.6,.8)*m);
                T density=1*unit_rho*scale_mass;
                Seed_Particles(box,0,0,density,particles_per_cell);
                int N_box_particles=particles.number-N_sphere_particles;
                ARRAY<int> foo(N_box_particles);
                for(int k=0;k<foo.m;k++) foo(k)=k+N_sphere_particles;
                Add_Fixed_Corotated(1*unit_p*scale_E,0.3,&foo);}

            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 11:{ // freefall circle, rising ground
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            ARRAY<int> foo(IDENTITY_ARRAY<>(particles.number));
            Add_Fixed_Corotated(1*unit_p*scale_E,0.3,&foo);
            Add_Gravity(m/(s*s)*TV(0,-9.8));
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1*m,false);
            MPM_COLLISION_OBJECT<TV>* bottom=dynamic_cast<MPM_COLLISION_OBJECT<TV>*>(collision_objects(3));
            bottom->func_frame=[this](T time){return FRAME<TV>(TV(0,(T).75-abs((T).15*time*scale_speed-(T).75)));};
            bottom->func_twist=[this](T time){return TWIST<TV>(TV(0,-sign((T).15*time*scale_speed-(T).75)*(T).15*scale_speed),typename TV::SPIN());};
        } break;
        case 12:{ // jello in a shrinking penalty box
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(TV(.3,.2)*m,TV(.7,.4)*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            ARRAY<int> foo(IDENTITY_ARRAY<>(particles.number));
            Add_Fixed_Corotated(20*unit_p*scale_E,0.4,&foo,false);
            Add_Gravity(m/(s*s)*TV(0,-9.8));
            Add_Walls(-1,COLLISION_TYPE::separate,1.9,.1*m,true);
            Add_Callbacks(true,"time-step",[this]()
                {
                    if(time>=10/24.0*s){
                        lagrangian_forces.Delete_Pointers_And_Clean_Memory();
                        this->solid_body_collection.deformable_body_collection.structures.Delete_Pointers_And_Clean_Memory();
                        this->output_structures_each_frame=true;
                        Add_Walls(-1,COLLISION_TYPE::separate,1.9,.1+(T)(time/s-10/24.0)*0.08*m,true);
                        Add_Gravity(m/(s*s)*TV(0,-9.8));}
                });
        } break;
        case 13:{ // surface tension test: fixed topology circle shell
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SEGMENTED_CURVE_2D<T>* sc=SEGMENTED_CURVE_2D<T>::Create();
            TV center(0.5*m,0.5*m);
            T radius=0.2*m;
            int N=80;
            T density=1*unit_rho*scale_mass;
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
            SEGMENTED_CURVE_2D<T>& new_sc=Seed_Lagrangian_Particles(*sc,[=](const TV& X){return TV(0.0,0)*(m/s);},0,density,true);
            SURFACE_TENSION_FORCE<TV>* stf=new SURFACE_TENSION_FORCE<TV>(new_sc,(T)0.01);
            Add_Force(*stf);
        } break;
        case 15:{ // colliding circles
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere0(TV(.25,.5)*m,.2*m);
            SPHERE<TV> sphere1(TV(.75,.5)*m,.2*m);
            T density=unit_rho*scale_mass;
            Seed_Particles(sphere0,[=](const TV& X){return TV(m/s,0);},0,density,particles_per_cell);
            Seed_Particles(sphere1,[=](const TV& X){return TV(-m/s,0);},0,density,particles_per_cell);
        } break;
        case 16:{ // oscillating circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            use_oldroyd=true;
            this->inv_Wi=(T)0;
            particles.Store_S(use_oldroyd);
            Seed_Particles(sphere,[=](const TV& X){return TV(0.1,0)*(m/s);},0,
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>()+1.5);
            particles.S.Fill(SYMMETRIC_MATRIX<T,2>()+sqr(1.5));
            OLDROYD_NEO_HOOKEAN<TV> *neo=new OLDROYD_NEO_HOOKEAN<TV>;
            neo->mu=38.462*unit_p*scale_E; // E=100, nu=0.3
            neo->lambda=57.692*unit_p*scale_E;
            MPM_OLDROYD_FINITE_ELEMENTS<TV> *fe=new MPM_OLDROYD_FINITE_ELEMENTS<TV>(force_helper,*neo,gather_scatter,0,this->inv_Wi,quad_F_coeff);
            Add_Force(*fe);
        } break;
        case 17:{ // spring test
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SEGMENTED_CURVE_2D<T>* sc=SEGMENTED_CURVE_2D<T>::Create();
            TV S(0.3*m,0.6*m);
            TV E(0.7*m,0.5*m);
            int N=20;
            TV D=(E-S)/(T)(N-1);
            T density=1*unit_rho*scale_mass;
            sc->particles.Add_Elements(N);
            sc->Update_Number_Nodes();
            for(int n=0;n<N-1;n++){
                TV X=S+D*n;
                sc->particles.X(n)=X;
                sc->mesh.elements.Append(TV_INT(n,n+1));}
            sc->particles.X(N-1)=E;
            SEGMENTED_CURVE_2D<T>& new_sc=Seed_Lagrangian_Particles(*sc,[=](const TV& X){return TV(0.0,0)*(m/s);},0,density,true);
            LINEAR_SPRINGS<TV>* stf=new LINEAR_SPRINGS<TV>(particles,new_sc.mesh);
            stf->Set_Restlength_From_Particles();
            stf->Set_Stiffness((T)10);
            stf->Set_Damping((T)0);
            Add_Force(*stf);
            for(int n=0;n<N;n++)
                for(int d=0;d<TV::m;d++)
                    particles.X(n)(d)+=random.Get_Uniform_Number(-0.2,0.2)*m;
        } break;
        case 18:{ // friction test
            T angle=use_foo_T1?foo_T1:0.1;
            T initial_velocity=use_foo_T2?foo_T2*(m/s):1*(m/s);
            T coefficient_of_friction=use_foo_T3?foo_T3:0.3;
            T g=9.8*m/(s*s);
            if(!no_regular_seeding) regular_seeding=true;
            Set_Grid(RANGE<TV>::Centered_Box()*m);
            RANGE<TV> box(TV(-.8,0)*m,TV(-.2,.2)*m);
            T density=unit_rho*scale_mass;
            Seed_Particles(box,[=](const TV& X){return TV(initial_velocity,0);},0,density,particles_per_cell);
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            FRAME<TV> frame(TV(),ROTATION<TV>::From_Angle(-angle));
            ORIENTED_BOX<TV> ground(RANGE<TV>(TV(-5,-5),TV(5,0))*m,frame);
            Add_Collision_Object(ground,COLLISION_TYPE::slip,coefficient_of_friction);
            for(int i=0;i<particles.X.m;i++){
                particles.X(i)=frame*particles.X(i);
                particles.V(i)=frame.r.Rotate(particles.V(i));}
            Add_Gravity(TV(0,-g));
            write_output_files.Append([=]()
                {
                    T c=sin(angle)*g;
                    T d=coefficient_of_friction*cos(angle)*g;
                    T v0=initial_velocity;
                    T v,x;
                    T t=time;
                    T st=FLT_MAX;
                    T a=c-sign(v0)*d;
                    if(v0*a<0) st=-v0/a;
                    if(t<st){
                        v=v0+a*t;
                        x=v0*t+.5*a*t*t;}
                    else{
                        x=.5*v0*st;
                        if(abs(c)>d){
                            a=c+sign(v0)*d;
                            t-=st;
                            v=a*t;
                            x+=.5*a*t*t;}
                        else v=0;}
                    Add_Debug_Particle(frame*(box.max_corner+TV(x,0)),VECTOR<T,3>(1,0,0));
                    Debug_Particle_Set_Attribute<TV>("V",frame.r.Rotate(TV(v,0)));
                });
        } break;
        case 20:{
            Set_Grid(RANGE<TV>(TV(0,0),TV(4,6))*m,TV_INT(2,3));
            int jet_freq=3;
            T init_vel=7*m;
            auto func=[this,jet_freq,init_vel](int frame)
                {
                    if (frame%jet_freq==0){
                        int old_m=particles.X.m;
                        T density=2*unit_rho*scale_mass;
                        RANGE<TV> box1(TV(.1,3.75)*m,TV(.3,3.77)*m), box2(TV(.1,3.5)*m,TV(.3,3.57));
                        Seed_Particles(box1,[=](const TV& X){return TV(init_vel,0);},0,density,particles_per_cell);
                        Seed_Particles(box2,[=](const TV& X){return TV(init_vel,0);},0,density,particles_per_cell);
                        box1+=RANGE<TV>(TV(3.6,0),TV(3.6,0))*m;
                        box2+=RANGE<TV>(TV(3.6,0),TV(3.6,0))*m;
                        Seed_Particles(box1,[=](const TV& X){return TV(-init_vel,0);},0,density,particles_per_cell);
                        Seed_Particles(box2,[=](const TV& X){return TV(-init_vel,0);},0,density,particles_per_cell);
                        ARRAY<int> new_particles(IDENTITY_ARRAY<>(particles.X.m-old_m)+old_m);
                        Set_Lame_On_Particles(1e3*unit_p*scale_E,.3,&new_particles);
                        particles.mu*=0;
                        particles.mu0*=0;
                    }
                };
            this->begin_frame.Append(func);
            Add_Collision_Object(RANGE<TV>(TV(0.,3.74),TV(.7,3.75))*m,COLLISION_TYPE::slip,0.);
            Add_Collision_Object(RANGE<TV>(TV(1.3,3.49),TV(2.,3.5))*m,COLLISION_TYPE::slip,0.);
            Add_Collision_Object(RANGE<TV>(TV(1.3,3.74),TV(2.,3.75))*m,COLLISION_TYPE::slip,0.);
            Add_Collision_Object(RANGE<TV>(TV(0.,3.49),TV(.7,3.5))*m,COLLISION_TYPE::slip,0.);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            Add_Walls(-1,COLLISION_TYPE::slip,.3,.1*m,false);
        } break;
        case 21:{ // circle with random initial velocities
            Set_Grid(RANGE<TV>(TV(-3,-3),TV(4,4))*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            random.Fill_Uniform(particles.V,-1*(m/s),1*(m/s));
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
        } break;
        case 22:{ // (fluid test) dam break reflection
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(grid.dX*2,TV(0.2,0.75)*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            if(extra_int.m && extra_int(0)) Add_Neo_Hookean(1e3*unit_p*scale_E,0.3);
            else Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            dilation_only=true;
            particles.mu*=0;
            particles.mu0*=0;
            Add_Gravity(m/(s*s)*TV(0,-1.8));
            Add_Walls(-1,COLLISION_TYPE::slip,0,0,false);
        } break;
        case 23:{ // (fluid test) dam break
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(grid.dX*2,TV(0.2,0.75)*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            dilation_only=true;
            particles.mu*=0;
            particles.mu0*=0;
            Add_Gravity(m/(s*s)*TV(0,-1.8));
            Add_Walls(-1,COLLISION_TYPE::slip,0,grid.dX.Min()*m,false);
        } break;
        case 24:{ // (fluid test) circle drop
            // one: ./mpm -kkt -scale_E 0
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.7)*m,.2*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 25:{ // (fluid test) pool of water w/ single particle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(grid.dX*4,TV(1*m-grid.dX(0)*4,0.25*m));
            T density=2*unit_rho*scale_mass;
            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Particle(TV(.5,.9),0,0,mass,volume);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            dilation_only=true;
            particles.mu*=0;
            particles.mu0*=0;
            Add_Walls(-1,COLLISION_TYPE::separate,.2,grid.dX(0)*4,false);
        } break;
        case 26:{ // Rayleigh Taylor
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            LOG::printf("FLAGS: %i %i\n",this->reflection_bc,this->reflection_bc_flags);
            int gap=this->reflection_bc?0:4;
            RANGE<TV> box(grid.dX*gap,TV(1*m-grid.dX(0)*gap,0.5*m));
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            ARRAY_VIEW<VECTOR<T,3> >* colors=particles.template Get_Array<VECTOR<T,3> >("color");
            for(int i=0;i<particles.X.m;i++)
                if(particles.X(i).y>box.Center().y)
                {
                    particles.mass(i)*=2;
                    (*colors)(i)=VECTOR<T,3>(1,0,0);
                }
            Add_Gravity(m/(s*s)*TV(0,-1.8));
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            dilation_only=true;
            particles.mu*=0;
            particles.mu0*=0;
            if(this->reflection_bc) this->reflection_bc_flags=-1;
            Add_Walls(-1,COLLISION_TYPE::separate,.2,grid.dX(0)*gap,false);
        } break;
        case 31:{ // Light fluid on heavy fluid
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(grid.dX*4,TV(1*m-grid.dX(0)*4,0.5*m));
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            ARRAY_VIEW<VECTOR<T,3> >* colors=particles.template Get_Array<VECTOR<T,3> >("color");
            for(int i=0;i<particles.X.m;i++)
                if(particles.X(i).y>box.Center().y)
                {
                    particles.mass(i)*=.5;
                    (*colors)(i)=VECTOR<T,3>(1,0,0);
                }
            Add_Gravity(m/(s*s)*TV(0,-1.8));
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            dilation_only=true;
            particles.mu*=0;
            particles.mu0*=0;
            Add_Walls(-1,COLLISION_TYPE::separate,.2,grid.dX(0)*4,false);
        } break;
        case 34:{ // Light fluid on heavy fluid
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(grid.dX*4,TV(1*m-grid.dX(0)*4,0.5*m));
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            ARRAY_VIEW<VECTOR<T,3> >* colors=particles.template Get_Array<VECTOR<T,3> >("color");
            
            int i=particles.X.m/3;
            particles.mass(i)*=.001;
            (*colors)(i)=VECTOR<T,3>(1,0,0);

            Add_Gravity(m/(s*s)*TV(0,-1.8));
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            dilation_only=true;
            particles.mu*=0;
            particles.mu0*=0;
            Add_Walls(-1,COLLISION_TYPE::separate,.2,grid.dX(0)*4,false);
        } break;
        case 98:{ // full box
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(grid.dX*(T)2,TV::All_Ones_Vector()*m-grid.dX*2);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 99:{ // single particle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            Add_Particle(TV(.8,.5),0,0,mass,volume);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 97:{ // single particle on collision box
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=1000*unit_rho*scale_mass;
            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            Add_Particle(TV(.5,.7),[](const TV& X){return TV(1,0);},0,mass,volume);
            RANGE<TV> box(TV(.1,.4)*m,TV(.9,.6)*m);
            Add_Collision_Object(box,COLLISION_TYPE::stick,friction);
            Add_Gravity(m/(s*s)*TV(0,-9.8));
        } break;
        case 27:{ // drop an oldroyd-b to a ground
            Set_Grid(RANGE<TV>(TV(-1,0),TV(1,1))*m,TV_INT(2,1));
            Add_Collision_Object(RANGE<TV>(TV(-5,-5),TV(5,.1))*m,COLLISION_TYPE::stick,0);
            SPHERE<TV> sphere(TV(.5,.5)*m,.2*m);
            T density=2*unit_rho*scale_mass;
            use_oldroyd=true;
            this->inv_Wi=(T)100;
            particles.Store_S(use_oldroyd);
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>()+1);particles.S.Fill(SYMMETRIC_MATRIX<T,2>()+sqr(1));
            LOG::cout<<particles.F<<std::endl<<std::endl;;
            LOG::cout<<particles.S<<std::endl;
            OLDROYD_NEO_HOOKEAN<TV> *neo=new OLDROYD_NEO_HOOKEAN<TV>;
            neo->mu=38.462*unit_p*scale_E; // E=100, nu=0.3
            neo->lambda=57.692*unit_p*scale_E;
            MPM_OLDROYD_FINITE_ELEMENTS<TV> *fe=new MPM_OLDROYD_FINITE_ELEMENTS<TV>(force_helper,*neo,gather_scatter,0,this->inv_Wi,quad_F_coeff);
            Add_Force(*fe);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 29:{ // drop an oldroyd-b to a ground SCA energy
            Set_Grid(RANGE<TV>(TV(-1,0),TV(1,1))*m,TV_INT(2,1));
            Add_Collision_Object(RANGE<TV>(TV(-5,-5),TV(5,.1))*m,COLLISION_TYPE::stick,0);
            SPHERE<TV> sphere(TV(.5,.5)*m,.2*m);
            T density=2*unit_rho*scale_mass;
            use_oldroyd=true;
            this->inv_Wi=(T)0;
            particles.Store_S(use_oldroyd);
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>()+1);particles.S.Fill(SYMMETRIC_MATRIX<T,2>()+sqr(1));
            VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV> *neo=new VOLUME_PRESERVING_OB_NEO_HOOKEAN<TV>;
            neo->mu=0;
            neo->lambda=57.692*unit_p*scale_E;
            MPM_OLDROYD_FINITE_ELEMENTS<TV> *fe=new MPM_OLDROYD_FINITE_ELEMENTS<TV>(force_helper,*neo,gather_scatter,0,this->inv_Wi,quad_F_coeff);
            Add_Force(*fe);
            Add_Force(*new MPM_VISCOSITY<TV>(force_helper,gather_scatter,0,0.001*unit_mu));
            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 28:{ // newton convergence problem: ./mpm 28 -affine -max_dt 1e-3 | grep converge
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            Add_Collision_Object(RANGE<TV>(TV(-5,-5),TV(5,.1))*m,COLLISION_TYPE::slip,0);
            Add_Collision_Object(RANGE<TV>(TV(-5,-5),TV(.1,5))*m,COLLISION_TYPE::slip,0);
            RANGE<TV> box(TV(.21,.21)*m,TV(.6,.5)*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
            Add_Fixed_Corotated(1.71*unit_p*scale_E,0.4);
        } break;
        case 30:{ // pinned rotating circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            VECTOR<T,1> angular_velocity(0.4*scale_speed/s);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,penalty_collisions_stiffness*kg/(s*s),
                penalty_damping_stiffness*kg/s);
            T pin_radius=sphere.radius*(T).8;
            for(int i=0;i<particles.X.m;i++)
                if((particles.X(i)-sphere.center).Magnitude_Squared()>sqr(pin_radius)){
                    TV dx=particles.X(i)-sphere.center;
                    pinning_force->Add_Target(i,
                        [=](T time){
                            ROTATION<TV> rot=ROTATION<TV>::From_Rotation_Vector(angular_velocity*time);
                            return rot.Rotate(dx)+sphere.center;});}
            Add_Force(*pinning_force);
            Add_Neo_Hookean(1e3*unit_p*scale_E,0.3);
        } break;
        case 32:{ // colliding
            Set_Grid(RANGE<TV>(TV(-1,-1),TV(2,2))*m);
            T d_small_1=1*unit_rho*scale_mass,d_small_2=1*unit_rho*scale_mass,d_large_1=d_small_1*.9,d_large_2=d_small_2*.9;
            SPHERE<TV> large1(TV(.2,.5)*m,.1*m);
            Seed_Particles(large1,[=](const TV& X){return TV(0.75,0)*(m/s);},0,d_large_1,particles_per_cell);
            SPHERE<TV> large2(TV(.6,.5)*m,.1*m);
            Seed_Particles(large2,0,0,d_large_2,particles_per_cell);
            for(int k=0;k<particles.number;k++)
                if((particles.X(k)-large1.center).Magnitude_Squared()<sqr(large1.radius*0.6)
                    || (particles.X(k)-large2.center).Magnitude_Squared()<sqr(large2.radius*0.6))
                    particles.deletion_list.Append(k);
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >("color");
            for(int i=0;i<particles.X.m;i++) (*color_attribute)(i)=VECTOR<T,3>(0,1,1);
            particles.Delete_Elements_On_Deletion_List();
            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            Add_Neo_Hookean(2*unit_p*scale_E,0.425,&mpm_particles);
            int old_pn=particles.number;
            SPHERE<TV> small1(large1.center,large1.radius*.6);
            Seed_Particles(small1,[=](const TV& X){return TV(0.75,0)*(m/s);},0,d_small_1,particles_per_cell);
            SPHERE<TV> small2(large2.center,large2.radius*.6);
            Seed_Particles(small2,0,0,d_small_2,particles_per_cell);
            ARRAY<int> mpm_particles2;for(int i=old_pn;i<particles.number;i++) mpm_particles2.Append(i);
            Add_Neo_Hookean(2000*unit_p*scale_E,0.3,&mpm_particles2);
            for(int i=old_pn;i<particles.X.m;i++) (*color_attribute)(i)=VECTOR<T,3>(1,1,1);
        } break;
        case 33:{ // sand box drop
            particles.Store_Fp(true);

            Set_Grid(RANGE<TV>::Unit_Box()*m);
            Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1))*m,COLLISION_TYPE::separate,10);

            T density=(T)1281*unit_rho*scale_mass;
            T E=5000*unit_p*scale_E,nu=.4;
            if(!use_theta_c) theta_c=0.01;
            if(!use_theta_s) theta_s=.00001;
            if(!use_hardening_factor) hardening_factor=80;
            if(!use_max_hardening) max_hardening=5;
            Add_Clamped_Plasticity(*new COROTATED_FIXED<T,TV::m>(E,nu),theta_c,theta_s,max_hardening,hardening_factor,0);
            RANGE<TV> box(TV(.45,.11)*m,TV(.55,.31)*m);
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-9.8));
        } break;
        case 35:{ // snow wedge
            // ./mpm 35 -flip 0.95 -max_dt .005 -cfl .1 -resolution 200
            particles.Store_Fp(true);

            Set_Grid(RANGE<TV>::Unit_Box()*m);
            ORIENTED_BOX<TV> wedge(RANGE<TV>(TV(),TV(0.2,0.2))*m,ROTATION<TV>::From_Angle(0.25*M_PI),TV(0.5,0.4-sqrt(2.0)*0.1));
            RANGE<TV> ground(TV(-1,0)*m,TV(2,0.1)*m);
            if(use_penalty_collisions){
                Add_Penalty_Collision_Object(wedge);
                Add_Penalty_Collision_Object(ground);}
            else{
                Add_Collision_Object(ground,COLLISION_TYPE::separate,1);
                Add_Collision_Object(new ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<TV> >(wedge),COLLISION_TYPE::separate,1);}

            T density=(T)2*unit_rho*scale_mass;
            int number_of_particles=20000;
            T E=40*unit_p*scale_E,nu=.2;
            if(!use_theta_c) theta_c=0.015;
            if(!use_theta_s) theta_s=.005;
            if(!use_hardening_factor) hardening_factor=7;
            Add_Clamped_Plasticity(*new COROTATED_FIXED<T,TV::m>(E,nu),theta_c,theta_s,max_hardening,hardening_factor,0);
            RANGE<TV> box(TV(.3,.7)*m,TV(.7,.9)*m);
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
            Add_Gravity(m/(s*s)*TV(0,-2));
        } break;
        case 36:{ // split
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,0,
                density,particles_per_cell);
            particles.F.Fill(MATRIX<T,2>(1,0,0,1));
            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            Add_Fixed_Corotated(unit_p*scale_E*10,0.3,&mpm_particles);
            PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,penalty_collisions_stiffness*kg/(s*s),penalty_damping_stiffness*kg/s);
            for(int i=0;i<particles.X.m;i++)
                if(particles.X(i)(0)<0.25*m || particles.X(i)(0)>90.75*m){TV x=particles.X(i);pinning_force->Add_Target(i,[=](T time){return x;});}
            Add_Force(*pinning_force);
            SEGMENTED_CURVE_2D<T>* sc=SEGMENTED_CURVE_2D<T>::Create();
            sc->particles.Add_Elements(particles.number);
            sc->Update_Number_Nodes();
            for(int p=0;p<particles.number;p++){
                sc->particles.X(p)=particles.X(p);
                for(int q=0;q<particles.number;q++){
                    TV L=particles.X(p),R=particles.X(q);
                    if(L(0)<=0.5*m && R(0)>0.5*m &&  L(1)>0.5*m && R(1)>0.5*m && (L-R).Magnitude_Squared()<sqr(grid.dX.Min()*2)){
                        sc->mesh.elements.Append(TV_INT(p,q));}}}
            LINEAR_SPRINGS<TV>* stf=new LINEAR_SPRINGS<TV>(particles,sc->mesh);
            ARRAY<T> r(sc->mesh.elements.m);r.Fill(grid.dX.Min()*2.1);
            stf->Set_Restlength(r);
            stf->Set_Stiffness((T)1); // TODO: units
            stf->Set_Damping((T)0.1); // TODO: units
            Add_Force(*stf);
            Add_Callbacks(true,"time-step",[this]()
                {
                    delete lagrangian_forces(lagrangian_forces.m-1);
                    lagrangian_forces.Remove_End();
                    SEGMENTED_CURVE_2D<T>* sc=SEGMENTED_CURVE_2D<T>::Create();
                    sc->particles.Add_Elements(particles.number);
                    sc->Update_Number_Nodes();
                    for(int p=0;p<particles.number;p++){
                        sc->particles.X(p)=particles.X(p);
                        for(int q=0;q<particles.number;q++){
                            TV L=particles.X(p),R=particles.X(q);
                            if(L(0)<=0.5*m && R(0)>0.5*m && L(1)>0.5*m && R(1)>0.5*m && (L-R).Magnitude_Squared()<sqr(grid.dX.Min()*2)){
                                sc->mesh.elements.Append(TV_INT(p,q));}}}
                    LINEAR_SPRINGS<TV>* stf=new LINEAR_SPRINGS<TV>(particles,sc->mesh);
                    ARRAY<T> r(sc->mesh.elements.m);r.Fill(grid.dX.Min()*2.1);
                    stf->Set_Restlength(r);
                    stf->Set_Stiffness((T)1); // TODO: units
                    stf->Set_Damping((T)0);
                    Add_Force(*stf);
                });
        } break;
        case 37:{ // sand box drop, better paramaters, with Hencky, usage: mpm 38 -resolution 100
            particles.Store_Fp(true);

            Set_Grid(RANGE<TV>::Unit_Box()*m);
            if(use_penalty_collisions)
                Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1))*m,0.9);
            else
                Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1))*m,COLLISION_TYPE::stick,0);

            T density=(T)1281*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.4;
            Add_St_Venant_Kirchhoff_Hencky_Strain(E,nu);
//            this->plasticity=new MPM_DRUCKER_PRAGER<TV>(friction_angle,cohesion);
            RANGE<TV> box(TV(.4,.1001)*m,TV(.45,.6001)*m);
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-9.81));
        } break;
        case 38:{ // sand box drop, wide
            particles.Store_Fp(true);

            Set_Grid(RANGE<TV>::Unit_Box()*m);
            if(use_penalty_collisions)
                Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1))*m,0.9);
            else
                Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1))*m,COLLISION_TYPE::stick,0);

            T density=(T)1281*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.4;
            Add_St_Venant_Kirchhoff_Hencky_Strain(E,nu);
//            this->plasticity=new MPM_DRUCKER_PRAGER<TV>(friction_angle,cohesion);
            T l0=0.05*m;
            T h0=l0*8;
            T gap=grid.dX(1)*0.01;
            RANGE<TV> box(TV(.5*m-l0,.1*m+gap),TV(.5*m+l0,.1*m+gap+h0));
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-9.81));
        } break;
        case 39:{ // DP on wedge
            particles.Store_Fp(true);

            Set_Grid(RANGE<TV>::Unit_Box()*m);
            ORIENTED_BOX<TV> wedge(RANGE<TV>(TV(),TV(0.2,0.2))*m,ROTATION<TV>::From_Angle(0.25*M_PI),TV(0.5,0.4-sqrt(2.0)*0.1)*m);
            RANGE<TV> ground(TV(-1,0)*m,TV(2,0.1)*m);
            if(use_penalty_collisions){
                Add_Penalty_Collision_Object(wedge);
                Add_Penalty_Collision_Object(ground);}
            else{
                Add_Collision_Object(ground,COLLISION_TYPE::separate,1);
                Add_Collision_Object(new ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<TV> >(wedge),COLLISION_TYPE::separate,1);}

            T density=(T)2*unit_rho*scale_mass;
            int number_of_particles=20000;
            T E=40*unit_p*scale_E,nu=.2;
            Add_St_Venant_Kirchhoff_Hencky_Strain(E,nu);
//            this->plasticity=new MPM_DRUCKER_PRAGER<TV>(friction_angle,cohesion);
            RANGE<TV> box(TV(.3,.7)*m,TV(.7,.9)*m);
            std::ifstream ifs("particles.dat");
            if(ifs.is_open()){
                std::string foobar;
                ifs>>foobar>>foobar>>foobar>>foobar;
                T x,y,mass,volume;
                while(ifs>>x>>y>>mass>>volume && !ifs.eof())
                    Add_Particle(TV(x,y)*m,0,0,mass*kg,volume*pow<TV::m>(m));
                PHYSBAM_ASSERT(particles.number==number_of_particles);}
            else{
                PHYSBAM_WARNING("Couldn't open 'particles.dat'. Falling back to using random particle positions.");
                Seed_Particles(box,0,0,density,number_of_particles*grid.dX.Product()/box.Size());}
            Add_Gravity(m/(s*s)*TV(0,-2));
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
        case 49:{ // sand column collapse (drucker prager v.s. snow-style)
            //  ./mpm 42 -use_exp_F -max_dt 7.5e-4 -scale_E 0.1 -resolution 100 -foo_T1 0(1)
            // fooT1=1:  use drucker prager
            // fooT1=0:  use snow-style plasticity
            particles.Store_Fp(true);
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            if(use_penalty_collisions) Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1))*m,0.9);
            else Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1))*m,COLLISION_TYPE::stick,0);
            T l0=0.05*m;
            T h0=l0*8;
            T gap=grid.dX(1)*0.01;
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            RANGE<TV> box(TV(.5*m-l0,.1*m+gap),TV(.5*m+l0,.1*m+gap+h0));
            if(foo_T1){
                if(!no_implicit_plasticity) use_implicit_plasticity=true;
                Seed_Particles(box,0,0,density,particles_per_cell);
                Add_Drucker_Prager_Case(E,nu,test_number-40);}
            else{
                if(!use_theta_c) theta_c=0.015;
                if(!use_theta_s) theta_s=.000001;
                if(!use_hardening_factor) hardening_factor=20;
                Add_Clamped_Plasticity(*new COROTATED_FIXED<T,TV::m>(E,nu),theta_c,theta_s,max_hardening,hardening_factor,0);
                Seed_Particles(box,0,0,density,particles_per_cell);}
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Drucker_Prager_Case(E,nu,test_number-40);
            Add_Gravity(m/(s*s)*TV(0,-9.81));
        } break;
        case 50:
        case 51:
        case 52:{ //lambda particles
            //usage:./mpm 51 -threads 8 -use_exp_F -max_dt 1e-3 -resolution 100 -scale_E 10 -fooT1 10 -fooT2 1000 -fooT3 4 -last_frame 20
            PHYSBAM_ASSERT(foo_T4<=particles_per_cell,"You can't have more water particles than total number of particles.");
            particles.Store_Fp(true);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            //particles.Store_Plastic_Deformation(true);
            Set_Grid(2*RANGE<TV>::Unit_Box()*m);
            if(use_penalty_collisions){
                Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(2.5,.1))*m,0.9);
                Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(.1,1.5))*m,0.9);}
            else{
                Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(2.5,.1))*m,COLLISION_TYPE::stick,0);
                Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(.1,1.5))*m,COLLISION_TYPE::stick,0);}

            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e5*unit_p*scale_E,nu=.3;
            T gap=grid.dX(1)*0.01;
            RANGE<TV> box(TV(.1*m+gap,.1*m+gap),TV(.3,.75)*m);
            //seed sand particles
            if(foo_T4<particles_per_cell) Seed_Particles(box,0,0,density,particles_per_cell-foo_T4);
            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager(E,nu,(T)35,&sand_particles);
            //seed lambda particles
            if(test_number==51){
                T El=500*foo_T1,nul=.1*foo_T3;
                Add_Lambda_Particles(&sand_particles,El,nul,foo_T2,true);}
            else if(test_number==52 && foo_T4!=0){
                T El=500*foo_T1,nul=.1*foo_T3;
                T lambdal=El*nul/((1+nul)*(1-2*nul));
                int ns=particles.X.m;
                Seed_Particles(box,0,0,foo_T2,foo_T4);
                ARRAY<int> lambda_particles(particles.X.m-ns);
                ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >("color");
                //fix volume and mass
                T volume=grid.dX.Product()/particles_per_cell;
                T mass_sand=density*volume;
                T mass_lambda=foo_T2*volume;
                for(int p=0;p<ns;p++){
                    particles.mass(p)=mass_sand;
                    particles.volume(p)=volume;}
                for(int k=0;k<lambda_particles.m;k++){
                    int p=ns+k;
                    lambda_particles(k)=p;
                    particles.mass(p)=mass_lambda;
                    particles.volume(p)=volume;
                    particles.mu(p)=(T)0;
                    particles.mu0(p)=(T)0;
                    particles.lambda(p)=lambdal;
                    particles.lambda0(p)=lambdal;
                    (*color_attribute)(p)=VECTOR<T,3>(0,0,1);}
                Add_Fixed_Corotated(El,nul,&lambda_particles,true);}

            Add_Gravity(m/(s*s)*TV(0,-9.81));
            if(extra_T.m>=2){
                MPM_VISCOSITY<TV>* visc=new MPM_VISCOSITY<TV>(force_helper,gather_scatter,0,0*unit_mu);
                visc->Use_Variable_Viscosity();
                for(int i=0;i<particles.X.m;i++)
                    visc->viscosity(i)=extra_T(particles.X(i).y>.4);
                Add_Force(*visc);}
        } break;
        case 53:
        case 56:
        case 57:{ // sandbox
            //  ./mpm 53 -threads 10 -use_exp_F -max_dt 7.5e-4 -scale_E 1 -resolution 200 -fooT1 10 -last_frame 10 -o sandbox_implicit
            particles.Store_Fp(true);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            if(use_penalty_collisions){
                LOG::cout<<"ERROR: NOT IMPLEMENTED FOR PENALTY."<<std::endl;;
                PHYSBAM_FATAL_ERROR();}
            else{
                // Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1))*m,COLLISION_TYPE::slip,0);
                RANGE<TV> ground(TV(-0.5,-1)*m,TV(2.5,.1)*m);
                RANGE<TV> leftwall(TV(-1,-1)*m,TV(.1,2.5)*m);
                RANGE<TV> rightwall(TV(0.9,-1)*m,TV(2.5,2.5)*m);
                RANGE<TV> top(TV(-0.5,.9)*m,TV(2.5,2)*m);
                Add_Collision_Object(
                    new IMPLICIT_OBJECT_UNION<TV>(
                        new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(ground),
                        new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(leftwall),
                        new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(rightwall),
                        new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(top)),
                    (test_number==53)?COLLISION_TYPE::slip:COLLISION_TYPE::stick,0);}

            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e5*unit_p*scale_E,nu=.3;
            T gap=grid.dX(1)*1.1;

            // SEED A REGULAR BOX
            // RANGE<TV> box(TV(.1*m+gap,.1*m+gap),TV(.9*m-gap,.3*m-gap));
            // Seed_Particles(box,0,0,density,particles_per_cell);

            // SEED FROM COLUMN COLLAPSES
            if(test_number==53){
                T ymax[5]={.2,.3,.24,.27,.4};
                T xmin[5]={.2,.33,.51,.63,.84};
                T xmax[5]={.3,.43,.59,.7,.88};
                for(int k=0;k<5;k++){
                    RANGE<TV> boxdune(TV(xmin[k]*m+gap,.1*m+gap),TV(xmax[k]*m-gap,ymax[k]*m-gap));
                    Seed_Particles(boxdune,0,0,density,particles_per_cell);}}
            else{
                RANGE<TV> box(TV(.1+gap,.1+gap)*m,TV(0.9-gap,.4)*m);
                Seed_Particles(box,0,0,density,particles_per_cell);
                //Add more collision object
                //foo_T2=stamp velocity
                //foo_T3=stamp original position
                if(!use_foo_T2) foo_T2=(T)0.6;
                if(!use_foo_T3) foo_T3=(T)0.5;
                TV min_corner((T)0.4,foo_T3);
                Add_Collision_Object(RANGE<TV>(min_corner*m,(min_corner+0.2)*m),COLLISION_TYPE::slip,0);}
            // SAND PROPERTIES
            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager(E,nu,(T)35,&sand_particles);
            Add_Gravity(m/(s*s)*TV(0,-9.81));
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >("color");
            for(int i=0;i<particles.X.m;i++) (*color_attribute)(i)=VECTOR<T,3>(.8,.7,.7);
            //Add water particles for case 57
            if(test_number==57){
                if(!use_foo_T4) foo_T4=(T)1;
                Add_Lambda_Particles(&sand_particles,E*foo_T4,nu,(T)1000*unit_rho,true,(T)0.3,(T)1);}

            // SEED an elastic dropping object
            // {int N_sand=particles.number;
            //     SPHERE<TV> sphere(TV(.5,.8)*m,0.05*m);
            //     T density=2900*unit_rho*foo_T1;
            //     Seed_Particles(sphere,0,0,density,particles_per_cell);
            //     int N_box_particles=particles.number-N_sand;
            //     ARRAY<int> foo(N_box_particles);
            //     for(int k=0;k<foo.m;k++) foo(k)=k+N_sand;
            //     Add_Fixed_Corotated(35.37e5*unit_p*scale_E,0.3,&foo);}
            if(test_number==56 || test_number==57)
                Add_Callbacks(true,"time-step",[this]()
                    {
                        T y=0;
                        T v=foo_T2*m/s;
                        T treshold=((T)0.2/v)*m;
                        T stop=2*treshold;
                        if(time<=stop){
                            lagrangian_forces.Delete_Pointers_And_Clean_Memory();
                            this->solid_body_collection.deformable_body_collection.structures.Delete_Pointers_And_Clean_Memory();
                            this->output_structures_each_frame=true;
                            delete collision_objects(collision_objects.m-1);
                            collision_objects.Pop();
                            if(time<=treshold) y=(foo_T3-time*v/s)*m;
                            else if(time<=stop) y=0.3;
                            else y=0.3+v*(time-2*treshold);
                            TV min_corner((T)0.4,y);
                            LOG::printf("time=%P\tmin_corner_y=%P\n",time,y);
                            Add_Collision_Object(RANGE<TV>(min_corner*m,(min_corner+0.2)*m),COLLISION_TYPE::slip,0);
                            Add_Gravity(m/(s*s)*TV(0,-9.8));}
                    });
        } break;
        case 54:
        case 58:
        case 61:{ // sand cup pull
            // ./mpm 54 -threads 10 -use_exp_F -max_dt 7.5e-4  -resolution 90 -last_frame 200
            particles.Store_Fp(true);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Set_Grid(RANGE<TV>(TV(-1,-0.5),TV(1,0.5))*m,TV_INT(2,1));
            RANGE<TV> cupbottom(TV(-0.3,-0.25),TV(0.3,-0.2));
            RANGE<TV> cupleft(TV(-0.3,-0.25),TV(-0.25,0.25));
            RANGE<TV> cupright(TV(0.25,-0.25),TV(0.3,0.25));
            Add_Collision_Object(
                new IMPLICIT_OBJECT_UNION<TV>(
                    new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(cupbottom),
                    new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(cupleft),
                    new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(cupright)),COLLISION_TYPE::separate,0);
            Add_Walls(-1,COLLISION_TYPE::stick,0,0.04,false);
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e5*unit_p*scale_E,nu=.3;
            ARRAY<RANGE<TV> > columns;
            columns.Append(RANGE<TV>(TV(-0.24,-0.2),TV(-0.2,0.1)));
            columns.Append(RANGE<TV>(TV(-0.15,-0.2),TV(-0.05,-0.1)));
            columns.Append(RANGE<TV>(TV(-0.02,-0.2),TV(0.1,0.02)));
            columns.Append(RANGE<TV>(TV(0.12,-0.2),TV(0.15,0.12)));
            columns.Append(RANGE<TV>(TV(0.18,-0.2),TV(0.23,-0.15)));
            for(int k=0;k<columns.m;k++) Seed_Particles(columns(k),0,0,density,particles_per_cell);
            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager(E,nu,(T)35,&sand_particles,false,test_number==61?sigma_Y:0);
            Add_Gravity(m/(s*s)*TV(0,-9.81));
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >("color");
            for(int i=0;i<particles.X.m;i++) (*color_attribute)(i)=VECTOR<T,3>(.8,.7,.7);
            //Add water particles for case 58
            if(test_number==58){
                if(!use_foo_T4) foo_T4=(T)1;
                Add_Lambda_Particles(&sand_particles,E*foo_T4,nu,(T)1000*unit_rho,true,(T)0.3,(T)1);}
            auto func=[this](int frame){if(frame==10) Add_Gravity(TV(0,20));};
            this->begin_frame.Append(func);
        } break;
        case 55:{ // Moving collision object
            particles.Store_Fp(true);

            Set_Grid(RANGE<TV>::Unit_Box()*m);
            if(use_penalty_collisions)
                PHYSBAM_NOT_IMPLEMENTED();
            Add_Collision_Object(RANGE<TV>(TV(0,.2),TV(1,.3))*m,COLLISION_TYPE::stick,0,
                [=](T time){return FRAME<TV>(TV(0,max(0.2/s*(time-.4),0.0))*m);},
                [=](T time){return TWIST<TV>(TV(0,0.2/s*(time-.4)>0?0.2*s:0)*m,typename TV::SPIN());});
            T density=(T)2200*unit_rho*scale_mass;
            T E=1e5*unit_p*scale_E,nu=.4;
            Add_Fixed_Corotated(E,nu);
            T l0=0.05*m;
            T h0=l0*8;
            T gap=grid.dX(1)*0.01;
            RANGE<TV> box(TV(.5*m-l0,.3*m+gap),TV(.5*m+l0,.3*m+gap+h0));
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-9.81));
        } break;
        case 59:{ // sand falling into a pile.
            particles.Store_Fp(true);
            Set_Grid(RANGE<TV>(TV(-2,0),TV(2,1))*m,TV_INT(4,1));
            RANGE<TV> ground(TV(-10,-10)*m,TV(10,.1)*m);
            if(use_penalty_collisions) Add_Penalty_Collision_Object(ground);
            else Add_Collision_Object(ground,COLLISION_TYPE::slip,0.3);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            int case_num=use_hardening_mast_case?hardening_mast_case:2;
            T E=1e4*unit_p*scale_E,nu=.3;
            Add_Drucker_Prager_Case(E,nu,case_num);
            TV spout(0,1*m);
            T density=(T)2200*unit_rho*scale_mass;
            T spout_width=.05*m;
            T pour_speed=.4*m/s;
            TV gravity=TV(0,-9.8*m/(s*s));
            Add_Source(spout,TV(0,-1),spout_width/2,pour_speed,gravity,density,E,nu,0,foo_T3);
            Add_Gravity(gravity);
        } break;
        case 60:{//cohesion sanity check
            particles.Store_Fp(true);

            Set_Grid(RANGE<TV>::Unit_Box()*m);
            if(use_penalty_collisions)
                Add_Penalty_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1))*m,0.9);
            else
                Add_Collision_Object(RANGE<TV>(TV(-0.5,-1),TV(1.5,.1))*m,COLLISION_TYPE::stick,0);
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            T l0=0.05*m;
            T h0=l0*8;
            T gap=grid.dX(1)*0.01;
            RANGE<TV> box(TV(.5*m-l0,.1*m+gap),TV(.5*m+l0,.1*m+gap+h0));
            Seed_Particles(box,0,0,density,particles_per_cell);
            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            if(!use_cohesion) sigma_Y=1;
            Add_Drucker_Prager(E,nu,(T)35,0,0,0,&sand_particles,false,sigma_Y);
            Add_Gravity(m/(s*s)*TV(0,-9.81));
        } break;
        case 62:{ // hourglass
            particles.Store_Fp(true);
            Set_Grid(RANGE<TV>(TV(-0.2,-0.45)*m,TV(0.2,0.45)*m),TV_INT(4,9),TV_INT(),4);
            LOG::cout<<"GRID DX: " <<grid.dX<<std::endl;
            IMPLICIT_OBJECT<TV>* hg=new ANALYTIC_IMPLICIT_OBJECT<HOURGLASS<TV> >(HOURGLASS<TV>(TV::Axis_Vector(1),TV(),(T).16,(T).0225,(T).8,(T).0225));
            IMPLICIT_OBJECT<TV>* inv=new IMPLICIT_OBJECT_INVERT<TV>(hg);
            if(use_penalty_collisions) Add_Penalty_Collision_Object(inv);
            else Add_Collision_Object(inv,COLLISION_TYPE::separate,.3);
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            RANGE<TV> fill_part=grid.domain;
            fill_part.min_corner.y=0;
            fill_part.max_corner.y=.1;
            ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> > io_fill_part(fill_part);
            IMPLICIT_OBJECT_INTERSECTION<TV> ioi(&io_fill_part,hg);
            ioi.owns_io.Fill(false);
            Seed_Particles(ioi,0,0,density,particles_per_cell);
            LOG::printf("added %i particles.",particles.X.m);
            Add_Drucker_Prager_Case(E,nu,2);
            Add_Gravity(m/(s*s)*TV(0,-9.81));
        } break;
        case 63:{
            particles.Store_Fp(true);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Set_Grid(RANGE<TV>(TV(-.05,-0.07),TV(0.05,.13))*m,TV_INT(1,2));
            RANGE<TV> ground(TV(-0.1,-2)*m,TV(2,-0.05)*m);
            RANGE<TV> leftwall(TV(-10,-10)*m,TV(-0.051,10)*m);
            RANGE<TV> rightwall(TV(0.049,-10)*m,TV(10,10)*m);
            Add_Collision_Object(ground,COLLISION_TYPE::slip,0.3);
            Add_Collision_Object(leftwall,COLLISION_TYPE::stick,0);
            Add_Collision_Object(rightwall,COLLISION_TYPE::stick,0);
            T cs=(T).25/20;
            RANGE<TV> cuptop(TV(-cs-0.01,0.025)*m,TV(cs+0.01,0.035)*m);
            RANGE<TV> cupleft(TV(-cs-0.01,-0.02)*m,TV(-cs,0.025)*m);
            RANGE<TV> cupright(TV(cs,-0.02)*m,TV(cs+0.01,0.025)*m);
            Add_Collision_Object(
                new IMPLICIT_OBJECT_UNION<TV>(
                    new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(cuptop),
                    new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(cupleft),
                    new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(cupright)),COLLISION_TYPE::separate,0);

            // voronoi
            TRIANGULATED_SURFACE<T>* strong=TRIANGULATED_SURFACE<T>::Create();
            Read_From_File(data_directory+"/../Private_Data/voronoi_strong_50.tri.gz",*strong);
            LOG::cout<<"Read mesh of strong #"<<strong->mesh.elements.m<<std::endl;
            LOG::cout<<"Read particles of strong #"<<strong->particles.number<<std::endl;
            for(int i=0;i<strong->particles.number;i++) strong->particles.X(i)/=20;
            strong->mesh.Initialize_Adjacent_Elements();
            strong->mesh.Initialize_Neighbor_Nodes();
            strong->mesh.Initialize_Incident_Elements();
            strong->Update_Bounding_Box();
            strong->Initialize_Hierarchy();
            strong->Update_Triangle_List();
            LOG::cout<<"Converting strong mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* strong_levelset=Initialize_Implicit_Surface(*strong,200);

            ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> > sphere(SPHERE<TV> (TV()*m,(.5/20)*m));
            RANDOM_NUMBERS<T> random;
            ARRAY<TV> X;
            POISSON_DISK<TV> poisson_disk(1);
            poisson_disk.Set_Distance_By_Volume(grid.dX.Product()/particles_per_cell);
            poisson_disk.Sample(random,sphere,X);

            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e5*unit_p*scale_E,nu=.3;
            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            for(int i=0;i<X.m;i++)
                if(strong_levelset->Extended_Phi(VECTOR<T,3>(X(i)(0),X(i)(1),0))<=(T)0)
                    Add_Particle(X(i),0,0,mass,volume);

            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager(E,nu,(T)35,&sand_particles,false,0);
            Add_Gravity(m/(s*s)*TV(0,9.81));
            //water
            if(!use_foo_T4) foo_T4=(T)1;
            Add_Lambda_Particles(&sand_particles,E*foo_T4,nu,(T)1000*unit_rho,true,(T)0.3,(T)1);
            int add_gravity_frame=restart?restart:60;
            auto func=[this,add_gravity_frame](int frame){if(frame==add_gravity_frame) Add_Gravity(TV(0,-20));};
            this->begin_frame.Append(func);
        } break;
        case 64:{
            particles.Store_Fp(true);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Set_Grid(RANGE<TV>(TV(-.04,-0.04),TV(0.04,.04))*m);
            RANGE<TV> ground(TV(-2,-2)*m,TV(2,-0.035)*m);
            RANGE<TV> leftwall(TV(-10,-10)*m,TV(-0.035,10)*m);
            RANGE<TV> rightwall(TV(0.035,-10)*m,TV(10,10)*m);
            Add_Collision_Object(
                new IMPLICIT_OBJECT_UNION<TV>(
                    new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(ground),
                    new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(leftwall),
                    new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(rightwall)),COLLISION_TYPE::slip,0.3);
            // voronoi
            TRIANGULATED_SURFACE<T>* strong=TRIANGULATED_SURFACE<T>::Create();
            Read_From_File(data_directory+"/../Private_Data/voronoi_strong_50.tri.gz",*strong);
            LOG::cout<<"Read mesh of strong #"<<strong->mesh.elements.m<<std::endl;
            LOG::cout<<"Read particles of strong #"<<strong->particles.number<<std::endl;
            for(int i=0;i<strong->particles.number;i++) strong->particles.X(i)/=20;
            strong->mesh.Initialize_Adjacent_Elements();
            strong->mesh.Initialize_Neighbor_Nodes();
            strong->mesh.Initialize_Incident_Elements();
            strong->Update_Bounding_Box();
            strong->Initialize_Hierarchy();
            strong->Update_Triangle_List();
            LOG::cout<<"Converting strong mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* strong_levelset=Initialize_Implicit_Surface(*strong,200);

            //full sphere
            TRIANGULATED_SURFACE<T>* full=TRIANGULATED_SURFACE<T>::Create();
            Read_From_File(data_directory+"/../Private_Data/voronoi_full_50.tri.gz",*full);
            LOG::cout<<"Read mesh of full #"<<full->mesh.elements.m<<std::endl;
            LOG::cout<<"Read particles of full #"<<full->particles.number<<std::endl;
            for(int i=0;i<full->particles.number;i++) full->particles.X(i)/=20;
            full->mesh.Initialize_Adjacent_Elements();
            full->mesh.Initialize_Neighbor_Nodes();
            full->mesh.Initialize_Incident_Elements();
            full->Update_Bounding_Box();
            full->Initialize_Hierarchy();
            full->Update_Triangle_List();
            LOG::cout<<"Converting the mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* full_levelset=Initialize_Implicit_Surface(*full,200);

            //sand
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e5*unit_p*scale_E,nu=.3;
            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            LOG::cout<<"Seeding particles..."<<std::endl;
            ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> > sphere(SPHERE<TV> (TV()*m,(.55/20)*m));
            RANDOM_NUMBERS<T> random;
            ARRAY<TV> X;
            POISSON_DISK<TV> poisson_disk(1);
            poisson_disk.Set_Distance_By_Volume(grid.dX.Product()/particles_per_cell);
            poisson_disk.Sample(random,sphere,X);
            for(int i=0;i<X.m;i++)
                if(full_levelset->Extended_Phi(VECTOR<T,3>(X(i)(0),X(i)(1),0))<=(T)0)
                    Add_Particle(X(i),[=](const TV& X){return TV(0,-0.7);},0,mass,volume);
            LOG::printf("Sand particles count=%P\n",particles.X.m);
            ARRAY<int> sand_particles(particles.X.m);
            for(int p=0;p<particles.X.m;p++) sand_particles(p)=p;
            Add_Drucker_Prager(E,nu,(T)35,&sand_particles,false,0);
            //water
            T porosity=0.3;
            T saturation_level=1;
            T water_density=(T)1000*unit_rho;
            T water_E=E*foo_T4;
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >("color");
            ARRAY<int> strong_lambda_particles;
            T volume_lambda=particles.volume(0)*porosity*saturation_level;
            T mass_lambda=water_density*volume_lambda;
            T lambda=water_E*nu/((1+nu)*(1-2*nu));
            for(int k=0;k<sand_particles.m;k++){
                if(strong_levelset->Extended_Phi(VECTOR<T,3>(particles.X(k)(0),particles.X(k)(1),0))<=0){
                    int p=particles.Add_Element();
                    particles.mass(p)=mass_lambda;
                    strong_lambda_particles.Append(p);
                    particles.lambda(p)=lambda;
                    particles.lambda0(p)=lambda;
                    (*color_attribute)(p)=VECTOR<T,3>(1,0,0);
                    (*color_attribute)(k)=VECTOR<T,3>(1,0,0);
                    int i=sand_particles(k);
                    particles.valid(p)=true;
                    particles.X(p)=particles.X(i);
                    particles.V(p)=particles.V(i);
                    particles.F(p)=particles.F(i);
                    if(particles.store_Fp) particles.Fp(p)=particles.Fp(i);
                    if(particles.store_B) particles.B(p)=particles.B(i);
                    if(particles.store_S) particles.S(p)=particles.S(i);
                    particles.volume(p)=volume_lambda;
                    particles.mu(p)=(T)0;
                    particles.mu0(p)=(T)0;}
                else{
                    (*color_attribute)(k)=VECTOR<T,3>(0,1,0);}}
            for(int p=0;p<particles.X.m;p++) particles.X(p)+=TV(0,-0.02);
            Add_Fixed_Corotated(water_E,nu,&strong_lambda_particles,true);
            Add_Gravity(m/(s*s)*TV(0,-9.81));
        } break;
        case 65:{
            //./mpm 65 -resolution 120 -max_dt 2e-6 -fooT4 10 -fooT2 1e-7 -cohesion 100 -no_implicit_plasticity -symplectic_euler -threads 8 -framerate 240 -last_frame 24
            particles.Store_Fp(true);
            particles.Store_Lame(true);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Set_Grid(RANGE<TV>(TV(-.06,-0.04),TV(0.06,.0))*m,TV_INT(3,1));
            RANGE<TV> ground(TV(-2,-2)*m,TV(2,-0.035)*m);
            RANGE<TV> leftwall(TV(-10,-10)*m,TV(-0.055,10)*m);
            RANGE<TV> rightwall(TV(0.055,-10)*m,TV(10,10)*m);
            Add_Collision_Object(
                new IMPLICIT_OBJECT_UNION<TV>(
                    new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(ground),
                    new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(leftwall),
                    new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(rightwall)),COLLISION_TYPE::slip,0.3);
            // voronoi
            TRIANGULATED_SURFACE<T>* strong=TRIANGULATED_SURFACE<T>::Create();
            Read_From_File(data_directory+"/../Private_Data/voronoi_strong_50.tri.gz",*strong);
            LOG::cout<<"Read mesh of strong #"<<strong->mesh.elements.m<<std::endl;
            LOG::cout<<"Read particles of strong #"<<strong->particles.number<<std::endl;
            for(int i=0;i<strong->particles.number;i++) strong->particles.X(i)/=20;
            strong->mesh.Initialize_Adjacent_Elements();
            strong->mesh.Initialize_Neighbor_Nodes();
            strong->mesh.Initialize_Incident_Elements();
            strong->Update_Bounding_Box();
            strong->Initialize_Hierarchy();
            strong->Update_Triangle_List();
            LOG::cout<<"Converting strong mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* strong_levelset=Initialize_Implicit_Surface(*strong,200);

            //full sphere
            TRIANGULATED_SURFACE<T>* full=TRIANGULATED_SURFACE<T>::Create();
            Read_From_File(data_directory+"/../Private_Data/voronoi_full_50.tri.gz",*full);
            LOG::cout<<"Read mesh of full #"<<full->mesh.elements.m<<std::endl;
            LOG::cout<<"Read particles of full #"<<full->particles.number<<std::endl;
            for(int i=0;i<full->particles.number;i++) full->particles.X(i)/=20;
            full->mesh.Initialize_Adjacent_Elements();
            full->mesh.Initialize_Neighbor_Nodes();
            full->mesh.Initialize_Incident_Elements();
            full->Update_Bounding_Box();
            full->Initialize_Hierarchy();
            full->Update_Triangle_List();
            LOG::cout<<"Converting the mesh to a level set..."<<std::endl;
            LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* full_levelset=Initialize_Implicit_Surface(*full,200);

            //sand
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e5*unit_p*scale_E,nu=.3;
            T mu=E/(2*(1+nu));
            T lambda=E*nu/((1+nu)*(1-2*nu));
            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            LOG::cout<<"Seeding particles..."<<std::endl;
            ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> > sphere(SPHERE<TV> (TV()*m,(.55/20)*m));
            RANDOM_NUMBERS<T> random;
            ARRAY<TV> X;
            POISSON_DISK<TV> poisson_disk(1);
            poisson_disk.Set_Distance_By_Volume(grid.dX.Product()/particles_per_cell);
            poisson_disk.Sample(random,sphere,X);
            for(int i=0;i<X.m;i++)
                if(full_levelset->Extended_Phi(VECTOR<T,3>(X(i)(0),X(i)(1),0))<=(T)0)
                    Add_Particle(X(i),[=](const TV& X){return TV(0,-2.7);},0,mass,volume);
            LOG::printf("Sand particles count=%P\n",particles.X.m);

            particles.mu.Fill(mu);
            particles.mu0.Fill(mu);
            particles.lambda.Fill(lambda);
            particles.lambda0.Fill(lambda);
            if(!use_foo_T2) foo_T2=1e-3;
            //if we are in the weak region, we use E that is scaled by foo_T2
            ARRAY<int> sand_particles;
            for(int k=0;k<particles.X.m;k++){
                sand_particles.Append(k);
                if(strong_levelset->Extended_Phi(VECTOR<T,3>(particles.X(k)(0),particles.X(k)(1),0))>0){
                    particles.mu(k)=mu*foo_T2;
                    particles.mu0(k)=mu*foo_T2;
                    particles.lambda(k)=lambda*foo_T2;
                    particles.lambda0(k)=lambda*foo_T2;}}
            Add_Drucker_Prager(0,0,(T)35,&sand_particles,false,sigma_Y);
            //water
            T porosity=0.3;
            T saturation_level=1;
            T water_density=(T)1000*unit_rho;
            if(!use_foo_T4) foo_T4=1e-3;
            T water_E=E*foo_T4;
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >("color");
            ARRAY<int> lambda_particles;
            T volume_lambda=particles.volume(0)*porosity*saturation_level;
            T mass_lambda=water_density*volume_lambda;
            T lambda_water_strong=water_E*nu/((1+nu)*(1-2*nu));
            T lambda_water_weak=foo_T2*water_E*nu/((1+nu)*(1-2*nu));
            for(int k=0;k<sand_particles.m;k++){
                bool inside=strong_levelset->Extended_Phi(VECTOR<T,3>(particles.X(k)(0),particles.X(k)(1),0))<=0;
                int p=particles.Add_Element();
                lambda_particles.Append(p);
                particles.valid(p)=true;
                particles.mass(p)=mass_lambda;
                int i=sand_particles(k);
                particles.X(p)=particles.X(i);
                particles.V(p)=particles.V(i);
                particles.F(p)=particles.F(i);
                if(particles.store_Fp) particles.Fp(p)=particles.Fp(i);
                if(particles.store_B) particles.B(p)=particles.B(i);
                if(particles.store_S) particles.S(p)=particles.S(i);
                particles.volume(p)=volume_lambda;
                particles.mu(p)=(T)0;
                particles.mu0(p)=(T)0;
                particles.lambda(p)=inside?lambda_water_strong:lambda_water_weak;
                particles.lambda0(p)=inside?lambda_water_strong:lambda_water_weak;
                (*color_attribute)(p)=inside?VECTOR<T,3>(1,0,0):VECTOR<T,3>(0,1,0);
                (*color_attribute)(k)=inside?VECTOR<T,3>(1,0,0):VECTOR<T,3>(0,1,0);}
            for(int p=0;p<particles.X.m;p++) particles.X(p)+=TV(0,-0.02);
            Add_Fixed_Corotated(water_E,nu,&lambda_particles,true);
            Add_Gravity(m/(s*s)*TV(0,-9.81));
        } break;

        case 66:{ // sand flow apic vs flip

            // ./mpm 66 -last_frame 36 -max_dt 2e-4 -scale_E 1000 -resolution 64 -fooT3 3 -particles_per_cell 20  -symplectic_euler -no_implicit_plasticity  -o sand_apic
            // ./mpm 66 -last_frame 36 -max_dt 2e-4 -scale_E 1000 -resolution 64 -fooT3 3 -particles_per_cell 20  -symplectic_euler -no_implicit_plasticity  -no_affine -flip 1 -o sand_flip

            particles.Store_Fp(true);
            Set_Grid(RANGE<TV>(TV(-2,0),TV(2,2))*m,TV_INT(4,2));
            RANGE<TV> ground(TV(-10,-10)*m,TV(10,.1)*m);

            if(use_penalty_collisions) Add_Penalty_Collision_Object(ground);
            else Add_Collision_Object(ground,COLLISION_TYPE::stick,0);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            int case_num=use_hardening_mast_case?hardening_mast_case:2;
            T E=1e4*unit_p*scale_E,nu=.3;
            Add_Drucker_Prager_Case(E,nu,case_num);
            TV spout(0,2*m);
            T density=(T)2200*unit_rho*scale_mass;
            T spout_width=.05*m;
            T pour_speed=.4*m/s;
            TV gravity=TV(0,-9.8*m/(s*s));
            Add_Source(spout,TV(0,-1),spout_width/2,pour_speed,gravity,density,E,nu,0,foo_T3);
            Add_Gravity(gravity);
        } break;

        case 67:{ // shake
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            Add_Walls(-1,COLLISION_TYPE::stick,0,.1*m,false);
            SEGMENTED_CURVE_2D<T>* sc=SEGMENTED_CURVE_2D<T>::Create();
            T density=1*unit_rho*scale_mass;
            int number_strainds=40;
            for(int z=0;z<number_strainds;z++){
                int newindex=sc->particles.X.m;
                TV S((0.4+z*0.2/number_strainds)*m,0.8*m);
                TV E((0.4+z*0.2/number_strainds)*m,0.3*m);
                int N=100;
                TV D=(E-S)/(T)(N-1);
                sc->particles.Add_Elements(N);
                sc->Update_Number_Nodes();
                for(int n=0;n<N-1;n++){
                    TV X=S+D*n;
                    sc->particles.X(newindex+n)=X;
                    sc->mesh.elements.Append(TV_INT(newindex+n,newindex+n+1));}
                sc->particles.X(newindex+N-1)=E;}
            SEGMENTED_CURVE_2D<T>& new_sc=Seed_Lagrangian_Particles(*sc,[=](const TV& X){return TV(0.0,0)*(m/s);},0,density,true);
            LINEAR_SPRINGS<TV>* stf=new LINEAR_SPRINGS<TV>(particles,new_sc.mesh);
            stf->Set_Restlength_From_Particles();
            stf->Set_Stiffness((T)50);
            stf->Set_Damping((T)0);
            Add_Force(*stf);

            ARRAY<int> mpm_particles(IDENTITY_ARRAY<>(particles.number));
            T E=10,nu=0.3;
            Add_Fixed_Corotated(unit_p*scale_E*E,nu,&mpm_particles,true);
            for(int i=0;i<particles.number;i++){
                particles.mu(i)=0;
                particles.mu0(i)=0;}

            PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,penalty_collisions_stiffness*kg/(s*s),
                penalty_damping_stiffness*kg/s);
            for(int i=0;i<particles.X.m;i++)
                if(particles.X(i)(1)>0.79){
                    TV dx=particles.X(i);
                    pinning_force->Add_Target(i,
                        [=](T time){
                            return dx+TV(0.2*sin(4*time),0);
                        });}
            Add_Force(*pinning_force);

            TV gravity=TV(0,-9.8*m/(s*s));
            Add_Gravity(gravity);
        } break;

        case 70:{ // Ringing test from original APIC paper
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            LOG::printf("DX %P %P\n",grid.dX,grid.one_over_dX);
            Add_Walls(10,COLLISION_TYPE::separate,.2,.1*m,false);
            Add_Gravity(TV(0,-2*m/(s*s)));
            max_dt=.005;
            cfl=0.1;

            // NOTE: mass is computed slightly differently.
            T mass_density=scale_mass*4*unit_rho;
            T E=scale_E*10*unit_p,nu=.3;
            LOG::printf("MASS DENSITY: %P\n",mass_density);

            RANGE<TV> range0(TV(.6,.1),TV(.8,.4));
            GRID<TV> grid0(TV_INT(TV(.625,.9375)*resolution)+1,range0);
            TRIANGULATED_AREA<T>& ta0=*new TRIANGULATED_AREA<T>;
            ta0.Initialize_Square_Mesh_And_Particles(grid0);
            TRIANGULATED_AREA<T>& new_ta0=Seed_Lagrangian_Particles(ta0,0,0,mass_density,true);
            Add_Fixed_Corotated(new_ta0,E,nu);

            T pinning_stiffness=1e5*unit_p;
            PINNING_FORCE<TV>* pinning_force=new PINNING_FORCE<TV>(particles,dt,pinning_stiffness,0);
            for(int i=0;i<grid0.counts.x;i++){
                TV X0=particles.X(i);
                pinning_force->Add_Target(i,[=](T time){return X0;});}
            for(int i=particles.number-grid0.counts.x;i<particles.number;i++){
                TV X0=particles.X(i);
                pinning_force->Add_Target(i,[=](T time){return X0;});}
            Add_Force(*pinning_force);

            RANGE<TV> range1(TV(.15,.2),TV(.25,.3));
            GRID<TV> grid1(TV_INT(TV(.3125,.3125)*resolution)+1,range1);
            TRIANGULATED_AREA<T>& ta1=*new TRIANGULATED_AREA<T>;
            ta1.Initialize_Square_Mesh_And_Particles(grid1);
            TRIANGULATED_AREA<T>& new_ta1=Seed_Lagrangian_Particles(ta1,[=](const TV& X){return TV(3*(m/s),0);},0,mass_density,true);
            Add_Fixed_Corotated(new_ta1,E,nu);
            LOG::printf("MASS: %P\n",particles.mass);

            for(int i=0;i<this->solid_body_collection.deformable_body_collection.structures.m;i++)
                this->solid_body_collection.deformable_body_collection.structures(i)->Update_Number_Nodes();
        } break;
        case 71:{ // (fluid test) dam break; Rabecca Brannon test
            Set_Grid(RANGE<TV>(TV(),TV(1,2)).Thickened(.25)*m,TV_INT(3,5),TV_INT(),2,100);
            if(this->reflection_bc_flags!=0)
                Set_Grid(RANGE<TV>(TV(),TV(1,2))*m,TV_INT(1,2),100);
            RANGE<TV> box(TV(.6,0)*m,TV(1,.4)*m);
            T density=1000*unit_rho;
            Seed_Particles(box,0,0,density,particles_per_cell);
            T stiffness=15e3*scale_E*unit_p,gamma=extra_T.m?extra_T(0):7;
            T tait_const=0.0894;
            ISOTROPIC_CONSTITUTIVE_MODEL<T,2>* cm=0;
            if(extra_int.m && extra_int(0)==1)
                cm=new TAIT_PRESSURE_FORCE<TV>(stiffness,tait_const);
            else cm=new QUASI_INCOMPRESSIBLE_FORCE<TV>(stiffness,gamma);
            DIAGONAL_MATRIX<T,2> dm;
            RANDOM_NUMBERS<T> rn;
            rn.Fill_Uniform(dm,.5,2);
            cm->Test(dm,0);
            Add_Force(*new MPM_FINITE_ELEMENTS<TV>(force_helper,*cm,gather_scatter,0));
            this->Add_Gravity2(m/(s*s)*TV(0,-9.81));
            if(this->reflection_bc_flags==0)
                Add_Walls(-1,COLLISION_TYPE::separate,0,.25*m,false);
            Add_Callbacks(false,"time-step",[=]()
                {
                    for(int i=0;i<particles.F.m;i++){
                        T J=particles.F(i).Determinant();
                        particles.F(i)=MATRIX<T,TV::m>()+pow<1,TV::m>(J);}
                    if(dt<1e-7) PHYSBAM_FATAL_ERROR("DT TOO SMALL");
                    for(int i=0;i<particles.V.m;i++)
                        if(particles.V(i).Magnitude()>50)
                            PHYSBAM_FATAL_ERROR("VELOCITY TOO BIG");
                });
        } break;
        case 72:{
            Set_Grid(RANGE<TV>(TV(),TV(1,2))*m,TV_INT(1,2));
            T density=2*unit_rho*scale_mass;
            T gravity=9.8*m/(s*s);
            T half_edge=.05*m;
            int seed_freq=4;
            int ob_per_frame=2;
            RANGE<TV> seed_box=grid.domain.Thickened(-.1*m);
            seed_box.min_corner.y=seed_box.max_corner.y-sqr(seed_freq*frame_dt)/2*gravity;
            seed_box=seed_box.Thickened(-sqrt((T)TV::m)*half_edge*m);
            seed_box.max_corner.x=seed_box.min_corner.x+seed_box.Edge_Lengths().x/ob_per_frame;
            auto func=[this,ob_per_frame,seed_freq,half_edge,density,seed_box](int frame)
                {
                    int num_part=particles.X.m;
                    if(frame%seed_freq==0 && frame<200){
                        for(int i=0;i<ob_per_frame;i++){
                            RANGE<TV> obj=RANGE<TV>::Centered_Box()*half_edge*m;
                            ROTATION<TV> rot;
                            random.Fill_Uniform(rot);
                            ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<TV> > ob(ORIENTED_BOX<TV>(obj,rot));
                            TV center;
                            random.Fill_Uniform(center,seed_box);
                            center.x+=i*seed_box.Edge_Lengths().x;
                            ob.analytic.corner+=center;
                            Seed_Particles(ob,0,0,density,particles_per_cell);}}
                    ARRAY<int> new_particles(IDENTITY_ARRAY<>(particles.X.m-num_part)+num_part);
                    Set_Lame_On_Particles(1e2*unit_p*scale_E,0.3,&new_particles);
                };
            this->begin_frame.Append(func);
            Add_Fixed_Corotated(1e2*unit_p*scale_E,0.3);
            Add_Gravity(TV(0,-gravity));
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1*m,false);
        } break;    

        case 68:{ // Source test
            Set_Grid(RANGE<TV>::Unit_Box()*m,TV_INT(1,1));
            T E=1e4*unit_p*scale_E,nu=.3;
            T density=(T)2200*unit_rho*scale_mass;
            TV gravity=TV(0,-9.8*m/(s*s));
            Add_Drucker_Prager_Case(E,nu,2);
            Add_Source(TV((T)0.5*m,0.9*m),TV(0,-1),.05*m/2,.4*m/s,gravity,density,E,nu,0,FLT_MAX);
            Add_Source(TV((T)0.9*m,(T)0.5*m),TV(-1,0),.1*m/2,4*m/s,gravity,density,E,nu,0,FLT_MAX);
            Add_Source(TV((T)0.2*m,(T)0.2*m),TV((T).4,(T).6),.05*m/2,3*m/s,gravity,density,E,nu,0,FLT_MAX);
            Add_Gravity(gravity);
            Add_Walls(-1,COLLISION_TYPE::separate,.3,.1*m,false);
        } break;    

        case 69:{ // Reflection bc test
            Set_Grid(RANGE<TV>(TV(),TV(4,1))*m,TV_INT(4,1));
            this->reflection_bc=true;
            this->use_full_reflection=true;
            T E=1e4*unit_p*scale_E,nu=.3;
            T density=(T)2*unit_rho*scale_mass;
            TV gravity=TV(0,-9.8*m/(s*s));
            T vel=1*m/s;
            RANGE<TV> box(TV(.25,0)*m,TV(1.25,.5)*m);
            VECTOR<TV,4> C;
            box.Corners(C);
            Seed_Particles(box,[=](const TV& X){return TV(vel,0);},0,density,particles_per_cell);
            Add_Fixed_Corotated(E,nu);
            Add_Gravity(gravity);
            Add_Walls(-1,COLLISION_TYPE::separate,.3,0,false);
            T a=this->reflection_bc_friction*-gravity.y;
            PHYSBAM_ASSERT(vel>=0);
            write_output_files.Append([this,vel,a,C]()
                {
                    TV v,x;
                    if(a*time<=vel)
                    {
                        v.x=vel-a*time;
                        x.x=vel*time-a*time*time/2;
                    }
                    else x.x=vel*vel/(a*2);
                    Add_Debug_Object<TV,2>({C(0)+x,C(1)+x},VECTOR<T,3>(1,0,0));
                    Add_Debug_Object<TV,2>({C(0)+x,C(2)+x},VECTOR<T,3>(1,0,0));
                    Add_Debug_Object<TV,2>({C(3)+x,C(1)+x},VECTOR<T,3>(1,0,0));
                    Add_Debug_Object<TV,2>({C(3)+x,C(2)+x},VECTOR<T,3>(1,0,0));
                });
        } break;    

        case 75:{ // Single Particel Stability Test
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            Add_Particle(TV(.52,.543),0,0,mass,volume);
            MATRIX<T,2> m;
            random.Fill_Uniform(m,-.1,.1);
            particles.F.Fill(m+1);
            T E=1e2*unit_p*scale_E;
            T nu=.3;
            Add_Neo_Hookean(E,nu);
        } break;

        case 76:{ // Reproduce the bowl instability problem
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            Add_Walls(-1,COLLISION_TYPE::separate,.1,.1*m,false);
            T density=4*unit_rho*scale_mass;
            Add_Collision_Object(SPHERE<TV>(TV(),.3*m),COLLISION_TYPE::separate,.1);
            Add_Collision_Object(SPHERE<TV>(TV(1,0)*m,.3*m),COLLISION_TYPE::separate,.1);
            SPHERE<TV> sphere1(TV(.1,.5)*m,.05*m);
            TRIANGULATED_AREA<T>* ta1=TESSELLATION::Generate_Triangles(sphere1,ceil(sphere1.radius/grid.dX.x*4));
            TRIANGULATED_AREA<T>& new_ta1=Seed_Lagrangian_Particles(*ta1,[=](const TV& X){return TV(foo_T1,0)*m;},
                0,density,true);
            Add_Fixed_Corotated(new_ta1,1e2*unit_p*scale_E,0.3);
            SPHERE<TV> sphere2(TV(.9,.5)*m,.05*m);
            TRIANGULATED_AREA<T>* ta2=TESSELLATION::Generate_Triangles(sphere2,ceil(sphere2.radius/grid.dX.x*4));
            TRIANGULATED_AREA<T>& new_ta2=Seed_Lagrangian_Particles(*ta2,[=](const TV& X){return TV(-foo_T1,0)*m;},
                0,density,true);
            Add_Fixed_Corotated(new_ta2,1e2*unit_p*scale_E,0.3);
            for(auto* s:this->solid_body_collection.deformable_body_collection.structures) s->Update_Number_Nodes();
            Add_Gravity(m/(s*s)*TV(0,-1.8));
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
