//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/PROJECTED_ARRAY.h>
#include <Core/Data_Structures/KD_TREE.h>
#include <Core/Matrices/FRAME.h>
#include <Core/Matrices/ROTATION.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Basic_Geometry/CONE.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/HOURGLASS.h>
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
#include <Geometry/Seeding/POISSON_DISK.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Rigids/Forces_And_Torques/RIGID_GRAVITY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Collisions_And_Interactions/PINNING_FORCE.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Constitutive_Models/MOONEY_RIVLIN_CURVATURE.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/OPENSUBDIV_SURFACE_CURVATURE_FORCE.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE_3D.h>
#include <Solids/Collisions/PENALTY_FORCE_COLLECTION.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_DRUCKER_PRAGER.h>
#include <Hybrid_Methods/Forces/MPM_OLDROYD_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_PLASTICITY_CLAMP.h>
#include <Hybrid_Methods/Forces/MPM_VISCOSITY.h>
#include <Hybrid_Methods/Forces/VOLUME_PRESERVING_OB_NEO_HOOKEAN.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include "STANDARD_TESTS_3D.h"
#include <fstream>
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
    if(!this->override_output_directory) viewer_dir.output_directory=LOG::sprintf("Test_3d_%i",test_number);
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
Write_Output_Files()
{
    if(write_output_files) write_output_files();
    BASE::Write_Output_Files();
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Read_Output_Files()
{
    if(read_output_files) read_output_files();
    BASE::Read_Output_Files();
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,3> >::
Initialize()
{
    this->asymmetric_system=true;
    switch(test_number)
    {
            // sticking: ./mpm_rb -3d 8 -float -scale_E .1 -rd_stiffness 1e-1 -max_dt .01 -T .1 -T .1 -T .2 -regular_seeding
            // sliding:
            // ./mpm_rb -3d 8 -float -scale_E .1 -rd_stiffness 1e-1 -max_dt .01 -T .1 -T .1 -T .2 -regular_seeding -rd_friction 0.1
        case 8:{ // MPM block and rigid circle, no collisions with boundary.
            T angle=extra_T(0);
            T vel=extra_T(1);
            T angle_t=extra_T(2);
            PHYSBAM_ASSERT(angle>=0 && angle<=pi/2 && angle_t>=0 && angle_t<=pi/2 && vel>=0);
            ROTATION<TV> Q(angle,TV(0,0,1)),Qn(angle_t,TV(0,1,0));
            ROTATION<TV> R=Q*Qn;
            TV n,t,t1;
            Q.Get_Rotated_Frame(t,n,t1);
            TV c(.5,.5,.5);

            Set_Grid(RANGE<TV>(TV(-1,0,0),TV(1,1,1))*m,TV_INT(2,1,1));

            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(R.Rotate(TV(0,0.1,0))+c,R));
            initial_state.twist.linear=-t*vel;
            GRID<TV> mattress_grid(TV_INT(10+1,10+1,10+1),RANGE<TV>(TV((T)-.1,(T)-.1,(T)-.1),TV((T).1,(T).1,(T).1))*m);
            T density=2*unit_rho*scale_mass;
            tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            particles.valid.Fill(true);
            particles.F.Fill(MATRIX<T,TV::m>()+1);
            if(particles.store_Fp) particles.Fp.Fill(MATRIX<T,TV::m>()+1);
            particles.volume.Fill(grid.domain.Size()/particles.number);
            
            for(int i=0;i<particles.number;i++){
                particles.V(i)=-t*vel;}
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >("color");
            color_attribute->Fill(TV(1,1,1));
            
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            TV g=m/(s*s)*TV(0,-1.8,0);
            Add_Gravity(g);
            Add_Collision_Object(Make_IO(PLANE<T>(n,c)));
            
            // Dump solution to viewer
            write_output_files=[=]()
            {
                T time=viewer_dir.frame_stack(0)*frame_dt;
                T mu=rd_penalty_friction;
                T acc=g.y*(cos(angle)*mu-sin(angle));
                if(vel+time*acc>0) // sliding
                {
                    T dist=(T).5*acc*sqr(time)+time*vel;
                    Add_Debug_Particle(c-t*dist,VECTOR<T,3>(0,1,0));
                }
                else // stopped
                {
                    T stop_time=-vel/acc;
                    T dist=(T).5*acc*sqr(stop_time)+stop_time*vel;
                    Add_Debug_Particle(c-t*dist,VECTOR<T,3>(1,1,0));
                }
                Add_Debug_Object(VECTOR<TV,2>(c-t,c+t),VECTOR<T,3>(1,0,0));
            };
        } break;

        case 9:{ // MPM block and static rigid plane
            T angle=extra_T(0);
            T vel=extra_T(1);
            T angle_t=extra_T(2);
            PHYSBAM_ASSERT(angle>=0 && angle<=pi/2 && angle_t>=0 && angle_t<=pi/2 && vel>=0);
            ROTATION<TV> Q(angle,TV(0,0,1)),Qn(angle_t,TV(0,1,0));
            ROTATION<TV> R=Q*Qn;
            TV n,t,t1;
            Q.Get_Rotated_Frame(t,n,t1);
            TV c(.5,.5,.5);

            Set_Grid(RANGE<TV>(TV(-1,0,0),TV(1,1,1))*m,TV_INT(2,1,1));

            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(R.Rotate(TV(0,0.1,0))+c,R));
            initial_state.twist.linear=-t*vel;
            GRID<TV> mattress_grid(TV_INT(10+1,10+1,10+1),RANGE<TV>(TV((T)-.1,(T)-.1,(T)-.1),TV((T).1,(T).1,(T).1))*m);
            T density=2*unit_rho*scale_mass;
            tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            particles.valid.Fill(true);
            particles.F.Fill(MATRIX<T,TV::m>()+1);
            if(particles.store_Fp) particles.Fp.Fill(MATRIX<T,TV::m>()+1);
            particles.volume.Fill(grid.domain.Size()/particles.number);
            
            for(int i=0;i<particles.number;i++){
                particles.V(i)=-t*vel;}
            ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >("color");
            color_attribute->Fill(TV(1,1,1));

            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            TV g=m/(s*s)*TV(0,-1.8,0);
            Add_Gravity(g);
            RIGID_BODY<TV>& rigid_body=tests.Add_Analytic_Box(TV(2,2,2));
            rigid_body.Frame().t=(c-R.Rotate(TV(0,1,0)))*m;
            rigid_body.Frame().r=R;
            rigid_body.is_static=true;

            // Dump solution to viewer
            write_output_files=[=]()
            {
                T time=viewer_dir.frame_stack(0)*frame_dt;
                T mu=rd_penalty_friction;
                T acc=g.y*(cos(angle)*mu-sin(angle));
                if(vel+time*acc>0) // sliding
                {
                    T dist=(T).5*acc*sqr(time)+time*vel;
                    Add_Debug_Particle(c-t*dist,VECTOR<T,3>(0,1,0));
                }
                else // stopped
                {
                    T stop_time=-vel/acc;
                    T dist=(T).5*acc*sqr(stop_time)+stop_time*vel;
                    Add_Debug_Particle(c-t*dist,VECTOR<T,3>(1,1,0));
                }
                Add_Debug_Object(VECTOR<TV,2>(c-t,c+t),VECTOR<T,3>(1,0,0));
            };
        } break;

        case 13:{ // Rigid-rigid version of inclined plane.
            T angle=extra_T(0);
            T vel=extra_T(1);
            T angle_t=extra_T(2);
            PHYSBAM_ASSERT(angle>=0 && angle<=pi/2 && angle_t>=0 && angle_t<=pi/2 && vel>=0);
            ROTATION<TV> Q(angle,TV(0,0,1)),Qn(angle_t,TV(0,1,0));
            ROTATION<TV> R=Q*Qn;
            TV n,t,t1;
            Q.Get_Rotated_Frame(t,n,t1);
            TV c(.5,.5,.5);

            Set_Grid(RANGE<TV>(TV(-2,0,0),TV(1,1,1))*m);
            RIGID_BODY<TV>& box=tests.Add_Analytic_Box(TV(.2,.2,.2),TV_INT()+1,(T)1);
            box.Frame().r=R;
            box.Frame().t=(R.Rotate(TV(0,.1,0))+c)*m;
            box.Twist().linear=-t*vel;
            TV g=m/(s*s)*TV(0,-1.8,0);
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
            RIGID_BODY<TV>& rigid_body=tests.Add_Analytic_Box(TV(2,2,2));
            rigid_body.Frame().t=(c-TV(0,1e-4,0)-Q.Rotate(TV(0,1,0)))*m;
            rigid_body.Frame().r=Q;
            rigid_body.is_static=true;

            // Dump solution to viewer
            write_output_files=[=]()
            {
                T time=viewer_dir.frame_stack(0)*frame_dt;
                T mu=rd_penalty_friction;
                T acc=g.y*(cos(angle)*mu-sin(angle));
                if(vel+time*acc>0) // sliding
                {
                    T dist=(T).5*acc*sqr(time)+time*vel;
                    Add_Debug_Particle(c-t*dist,VECTOR<T,3>(0,1,0));
                }
                else // stopped
                {
                    T stop_time=-vel/acc;
                    T dist=(T).5*acc*sqr(stop_time)+stop_time*vel;
                    Add_Debug_Particle(c-t*dist,VECTOR<T,3>(1,1,0));
                }
                Add_Debug_Object(VECTOR<TV,2>(c-t,c+t),VECTOR<T,3>(1,0,0));
            };
        } break;

            // Diff test for IO-MPM penalty force
            // ./mpm_rb -3d 201 -double -rd_stiffness 1e2 -test_diff
        case 201:{
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            Add_Collision_Object(Make_IO(PLANE<T>(TV(0,1,0),TV(0.5,0.5,0.5))));

            write_output_files=[]()
            {
                Add_Debug_Object(VECTOR<TV,3>(
                    TV(0,0.5,0),TV(1,0.5,0),TV(1,0.5,1)),VECTOR<T,3>(1,0,0));
                Add_Debug_Object(VECTOR<TV,3>(
                    TV(1,0.5,1),TV(0,0.5,1),TV(0,0.5,0)),VECTOR<T,3>(1,0,0));
            };

            T density=2*unit_rho*scale_mass;
            T volume=grid.dX.Product();
            T mass=density*volume;
            Add_Particle(TV(0.5,0.6,0.5),0,0,mass,volume);
            particles.V(0)=TV(0.1,-0.1,0.1);
        } break;

            // Diff test for Rigid-MPM penalty force
            // ./mpm_rb -3d 202 -double -rd_stiffness 1e1 -test_diff
        case 202:{
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            RIGID_BODY<TV>& cy=tests.Add_Analytic_Cylinder((T)1,(T)0.5,16,32,density);
            cy.Frame().t=TV(0.5,0,0.5);

            T volume=grid.dX.Product();
            T mass=density*volume;
            Add_Particle(TV(0.5,0.6,0.5),0,0,mass,volume);
            particles.V(0)=TV(0.1,-0.1,0.1);
        } break;

            //  ./mpm_rb -3d 101 -double -rd_stiffness 1e-1 -max_dt .01 -T 0.6 -rd_friction 0.5 -last_frame 200
        case 101:{
            T angle=extra_T(0);
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RIGID_BODY<TV>& ground=tests.Add_Rigid_Body("ground",(T)1,(T)0);
            ground.Frame().t=TV(0,0.1,0);
            ground.is_static=true;
            T density=10;
            RIGID_BODY<TV>& cy1=tests.Add_Analytic_Cylinder((T)0.2,(T)0.01,16,32,density);
            cy1.Frame().t=TV(0,0.11,0)*m;
            RIGID_BODY<TV>& cy2=tests.Add_Analytic_Cylinder((T)0.2,(T)0.01,16,32,density);
            cy2.Frame().r=ROTATION<TV>(angle,TV(0,1,0));
            cy2.Frame().t=(TV(0,0.13,0)+cy2.Frame().r.Rotate(TV(0,0,0.05)))*m;
            TV g=m/(s*s)*TV(0,-1.8,0);
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
        } break;

            // SIGGRAPH test 7, Wedging test, Rigid-MPM
            // sticking: ./mpm_rb -3d 71 -double -scale_E .1 -rd_stiffness 1e2 -max_dt .01 -regular_seeding -rd_friction 0.25
            // slipping: ./mpm_rb -3d 71 -double -scale_E .1 -rd_stiffness 1e2 -max_dt .01 -regular_seeding -rd_friction 0.1
        case 71:{
            Set_Grid(RANGE<TV>::Centered_Box()*2*m);
            T density=2*unit_rho*scale_mass;
            RIGID_BODY<TV>& lw=tests.Add_Analytic_Box(TV(0.4,8,8));
            lw.Frame().t=TV(-0.2,0.5,0.5);
            lw.is_static=true;
            Add_Collision_Object(Make_IO(PLANE<T>(TV(-1,0,0),TV(1,0.5,0.5))));

            TV g=m/(s*s)*TV(0,-1.8,0);
            RIGID_BODY<TV>& cube=tests.Add_Analytic_Box(TV(0.4,0.4,0.4),TV_INT()+1,density);
            cube.Frame().t=TV(0.2,0.5,0.5)*m;
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);

            SPHERE<TV> sphere(TV(.7,.5,0.5)*m,.4*m);
            TRIANGULATED_SURFACE<T>* ts=TESSELLATION::Generate_Triangles(sphere,4);
            Seed_Particles_Surface(*ts,*Make_IO(sphere),0,0,density,particles_per_cell);
            for(int i=0;i<solid_body_collection.deformable_body_collection.structures.m;i++)
                solid_body_collection.deformable_body_collection.structures(i)->Update_Number_Nodes();
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            Add_Gravity(g);
        } break;

            // SIGGRAPH test 7, Wedging test, Rigid-Rigid
            // sticking: ./mpm_rb -3d 72 -double -rd_stiffness 1e1 -max_dt .01 -rd_friction 0.375
            // slipping: ./mpm_rb -3d 72 -double -rd_stiffness 1e1 -max_dt .01 -rd_friction 0.3
        case 72:{
            Set_Grid(RANGE<TV>::Centered_Box()*2*m);
            RIGID_BODY<TV>& lw=tests.Add_Analytic_Box(TV(0.4,8,8));
            lw.Frame().t=TV(-0.2,0.5,0.5);
            lw.is_static=true;
            RIGID_BODY<TV>& rw=tests.Add_Analytic_Box(TV(0.4,8,8));
            rw.Frame().t=TV(1.2,0.5,0.5);
            rw.is_static=true;

            TV g=m/(s*s)*TV(0,-1.8,0);
            T density=1;
            RIGID_BODY<TV>& lcube=tests.Add_Analytic_Box(TV(0.52,0.5,0.5),TV_INT()+5,density);
            lcube.Frame().t=TV(0.25,0.5,0.5)*m;
            RIGID_BODY<TV>& rcube=tests.Add_Analytic_Box(TV(0.52,0.5,0.5),TV_INT()+5,density);
            rcube.Frame().t=TV(0.75,0.5,0.5)*m;

            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
        } break;

            // SIGGRAPH test 12, rolling rigid ball
            // rolling: ./mpm_rb -3d 12 -double -rd_stiffness 1e1 -max_dt .01 -T .2 -rd_friction 0.3
            // slipping: ./mpm_rb -3d 12 -double -rd_stiffness 1e1 -max_dt .01 -T .2 -rd_friction 0.1
        case 12:{
            T angle=extra_T(0);
            PHYSBAM_ASSERT(angle>=0 && angle<=pi/2);
            ROTATION<TV> Q(angle,TV(0,0,1));
            TV n,t,t1;
            Q.Get_Rotated_Frame(t,n,t1);
            TV c(.8,.5,.5);

            RIGID_BODY<TV>& sphere=tests.Add_Rigid_Body("sphere",(T)0.2,(T)0);
            sphere.Frame().t=c+0.2*n;

            Set_Grid(RANGE<TV>(TV(-2,0,0),TV(1,1,1))*m);
            TV g=m/(s*s)*TV(0,-9.8,0);
            RIGID_BODY<TV>& rigid_body=tests.Add_Analytic_Box(TV(4,2,2));
            rigid_body.Frame().t=(c-Q.Rotate(TV(0,1,0)))*m;
            rigid_body.Frame().r=Q;
            rigid_body.is_static=true;
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
        } break;

        case 60:{ // sand on obj
            particles.Store_Fp(true);
            Set_Grid(RANGE<TV>(TV(),TV(3,1,3))*m,TV_INT(3,1,3));

            RIGID_BODY<TV>& bw=tests.Add_Analytic_Box(TV(3,0.1,3));
            bw.is_static=true;
            bw.Frame().t=TV(1.5,0.05,1.5);

            use_colored_sand=true;
            T E=35.37e4*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Add_Drucker_Prager(E,nu,foo_T2);
            TV spout(1.5,0.8,1.5);
            T density=(T)2200*unit_rho*scale_mass;
            T spout_width=.1*m;
            T pour_speed=.2*m/s;
            TV gravity=TV(0,-9.8*m/(s*s),0);
            Set_Lame_On_Particles(E,nu);
            Add_Source(spout,TV(0,-1,0),spout_width/2,pour_speed,gravity,density,E,nu,0,FLT_MAX);
            Add_Gravity(gravity);

            RIGID_BODY<TV>& sphere=tests.Add_Analytic_Sphere(0.2,density*0.1);
            sphere.Frame().t=TV(1.35,0.3,1.35);
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,gravity);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
        } break;

        case 61:{ // goo on obj
            particles.Store_Fp(true);
            particles.Store_Lame(true);
            particles.Store_Lame0(true);
            Set_Grid(RANGE<TV>(TV(),TV(3,1,3))*m,TV_INT(3,1,3));

            RIGID_BODY<TV>& bw=tests.Add_Analytic_Box(TV(3,0.1,3));
            bw.is_static=true;
            bw.Frame().t=TV(1.5,0.05,1.5);

            T E=1e4*unit_p*scale_E,nu=.3;
            TV spout(1.5,0.8,1.5);
            T density=(T)2200*unit_rho*scale_mass;
            T spout_width=.1*m;
            T pour_speed=.2*m/s;
            TV gravity=TV(0,-9.8*m/(s*s),0);
            Set_Lame_On_Particles(E,nu);
            Add_Source(spout,TV(0,-1,0),spout_width/2,pour_speed,gravity,density,E,nu,0,FLT_MAX);
            Add_Gravity(gravity);

            if(!use_theta_c) theta_c=0.015;
            if(!use_theta_s) theta_s=.000001;
            if(!use_hardening_factor) hardening_factor=20;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Add_Clamped_Plasticity(*new COROTATED_FIXED<T,TV::m>(E,nu),theta_c,theta_s,0,hardening_factor,0);
            Set_Lame_On_Particles(E,nu);
            Add_Gravity(gravity);

            RIGID_BODY<TV>& sphere=tests.Add_Analytic_Sphere(0.2,density*0.1);
            sphere.Frame().t=TV(1.35,0.3,1.35);
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,gravity);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
        } break;

        case 21:{
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            T gap=0.005;
            T r=0.15;
            T R=r+gap;
            T k=0.8;
            SPHERE<TV> el(TV(0.5-R,r,0.5+tan(pi/6)*R),r);
            TRIANGULATED_SURFACE<T>* ts=TESSELLATION::Generate_Triangles(el,4);
            Seed_Particles_Surface(*ts,*Make_IO(el),0,0,density,particles_per_cell);
            ts->particles.X+=TV(2*R,0,0);
            el.center+=TV(2*R,0,0);
            Seed_Particles_Surface(*ts,*Make_IO(el),0,0,density,particles_per_cell);
            ts->particles.X+=TV(-R,0,-sin(pi/3)*2*R);
            el.center+=TV(-R,0,-sin(pi/3)*2*R);
            Seed_Particles_Surface(*ts,*Make_IO(el),0,0,density,particles_per_cell);

            for(int i=0;i<solid_body_collection.deformable_body_collection.structures.m;i++)
                solid_body_collection.deformable_body_collection.structures(i)->Update_Number_Nodes();
            
            particles.X.template Project<T,&TV::y>()*=k;
            particles.X.template Project<T,&TV::y>()+=0.1;
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            TV g=m/(s*s)*TV(0,-1.8,0);
            Add_Gravity(g);
            RIGID_BODY<TV>& sphere=tests.Add_Analytic_Sphere(r,15*density);
            sphere.Frame().t=TV(0.5,r*k+sqrt(6.0)/3*2*r+0.1-0.03,0.5);
            RIGID_BODY<TV>& ground=tests.Add_Analytic_Box(TV(1,0.1,1));
            ground.Frame().t=TV(0.5,0.05,0.5);
            ground.is_static=true;
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
        } break;

        case 720:{
            Set_Grid(RANGE<TV>::Centered_Box()*10*m);
            T density=1000;
            tests.Add_Analytic_Box(TV()+2,TV_INT()+resolution,density).Frame().t.y=7.04*m;
            tests.Add_Analytic_Box(TV()+2,TV_INT()+resolution,density).Frame().t.y=9.05*m;
            solid_body_collection.rigid_body_collection.rigid_body_particles.twist(0).linear.y=1;
            solid_body_collection.rigid_body_collection.rigid_body_particles.twist(1).linear.y=-1;
            break;}

        case 140:{
            particles.Store_Fp(true);
            Set_Grid(RANGE<TV>(TV(-1,-1,-1),TV(1,3,1))*m);

            RANDOM_NUMBERS<T> rng(seed);
            T density=(T)2200*unit_rho*scale_mass;
            TV g=m/(s*s)*TV(0,1.8,0);
            RIGID_BODY<TV>& bowl=tests.Add_Analytic_Bowl(TV(),TV(0,1,0),(T)0,(T)1,(T)0.05);
            bowl.is_static=true;

            RANGE<TV> box(TV(-0.6,-2.6,-0.6),TV(0.6,0,0.6));
            RANGE<TV> R(TV(0,0,0),TV(pi,pi,pi));
            POISSON_DISK<TV> poisson_disk(1);
            ARRAY<TV> X;
            T gap=0.45;
            poisson_disk.Set_Distance_By_Volume(cube(gap));
            poisson_disk.Sample(rng,box,X);
            for(int i=0;i<X.m/2;i++){
                RIGID_BODY<TV>& ring=tests.Add_Rigid_Body("Rings_Test/ring_revolve",(T)0.05,(T)0);
                ring.Set_Mass(ring.Volume()*density);
                TV rotation;
                rng.Fill_Uniform(rotation,R);
                ring.Frame().r=ROTATION<TV>::From_Euler_Angles(rotation.x,rotation.y,rotation.z);
                ring.Frame().t=X(i);}

            for(int i=X.m/2;i<X.m;i++){
                RIGID_BODY<TV>& torus=tests.Add_Analytic_Torus((T)0.05,(T)0.1,16,32,density);
                TV rotation;
                rng.Fill_Uniform(rotation,R);
                torus.Frame().r=ROTATION<TV>::From_Euler_Angles(rotation.x,rotation.y,rotation.z);
                torus.Frame().t=X(i);}

            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);

            begin_frame=[=](int frame)
            {
                if(frame!=45) return;
                RANGE<TV> sand_cube(TV(-0.4,-0.3,-0.4),TV(0.4,0,0.4));
                T E=35.37e6*unit_p*scale_E,nu=.3;
                if(!use_theta_c) theta_c=0.015;
                if(!use_theta_s) theta_s=.000001;
                if(!use_hardening_factor) hardening_factor=20;
                Add_Clamped_Plasticity(*new COROTATED_FIXED<T,TV::m>(E,nu),theta_c,theta_s,max_hardening,hardening_factor,0);
                Seed_Particles(sand_cube,0,0,density,particles_per_cell);
                ARRAY_VIEW<VECTOR<T,3> >* colors=particles.template Get_Array<VECTOR<T,3> >("color");
                for(int p=0;p<particles.number;p++) (*colors)(p)=Sand_Color();
                Add_Drucker_Prager_Case(E,nu,2);
                Add_Gravity(g);
            };
            break;}

        case 150:{
            particles.Store_Fp(true);
            Set_Grid(RANGE<TV>(TV(-1.2,-1.2,-1.2),TV(1.2,1.2,1.2))*m);

            RANDOM_NUMBERS<T> rng(seed);
            T density=(T)2200*unit_rho*scale_mass;
            TV g=m/(s*s)*TV(0,-1.8,0);
            RIGID_BODY<TV>& bowl=tests.Add_Analytic_Bowl(TV(),TV(0,1,0),(T)0,(T)1,(T)0.05);
            bowl.is_static=true;

            RANGE<TV> box(TV(-0.6,0,-0.6),TV(0.6,5.6,0.6));
            RANGE<TV> R(TV(0,0,0),TV(pi,pi,pi));
            POISSON_DISK<TV> poisson_disk(1);
            ARRAY<TV> X;
            T gap=0.5;
            poisson_disk.Set_Distance_By_Volume(cube(gap));
            poisson_disk.Sample(rng,box,X);
            rng.Random_Shuffle(X);
            auto intersects=[this](const RIGID_BODY<TV>& a,const RIGID_BODY<TV>& b)
            {
//                if(!a.Bounding_Boxes_Intersect(b)) return false;
                for(int i=0;i<a.simplicial_object->particles.X.m;i++)
                    if(b.Implicit_Geometry_Lazy_Inside(a.World_Space_Point(a.simplicial_object->particles.X(i)),.02))
                        return true;
                for(int i=0;i<b.simplicial_object->particles.X.m;i++)
                    if(a.Implicit_Geometry_Lazy_Inside(b.World_Space_Point(b.simplicial_object->particles.X(i)),.02))
                        return true;
                return false;
            };
            auto intersects_any=[this,intersects](const RIGID_BODY<TV>& a)
            {
                for(int i=0;i<solid_body_collection.rigid_body_collection.rigid_body_particles.number;i++)
                    if(i!=a.particle_index)
                        if(intersects(a,solid_body_collection.rigid_body_collection.Rigid_Body(i)))
                            return true;
                return false;
            };
            
            for(int i=0;i<X.m/2;i++){
                RIGID_BODY<TV>& ring=tests.Add_Rigid_Body("Rings_Test/ring_revolve",(T)0.05*1.5,(T)0);
                ring.Set_Mass(ring.Volume()*density);
                TV rotation;
                for(int j=0;j<20;j++){
                    LOG::printf("try %i %i\n",i,j);
                    rng.Fill_Uniform(rotation,R);
                    ring.Frame().r=ROTATION<TV>::From_Euler_Angles(rotation.x,rotation.y,rotation.z);
                    ring.Frame().t=X(i);
                    if(!intersects_any(ring)) break;
                    PHYSBAM_ASSERT(j!=19);}}

            for(int i=X.m/2;i<X.m;i++){
                RIGID_BODY<TV>& torus=tests.Add_Analytic_Torus((T)0.05*1.5,(T)0.1*1.5,16,32,density);
                for(int j=0;j<20;j++){
                    LOG::printf("try %i %i\n",i,j);
                    TV rotation;
                    rng.Fill_Uniform(rotation,R);
                    torus.Frame().r=ROTATION<TV>::From_Euler_Angles(rotation.x,rotation.y,rotation.z);
                    torus.Frame().t=X(i);
                    if(!intersects_any(torus)) break;
                    PHYSBAM_ASSERT(j!=19);}}

            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);

            use_colored_sand=true;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Add_Drucker_Prager(E,nu,foo_T2);
            TV spout(0,1,0);
            T spout_width=.2*m;
            T pour_speed=.3*m/s;
            T source_start=4;
            T source_end=100;
            Set_Lame_On_Particles(E,nu);
            Add_Source(spout,TV(0,-1,0),spout_width/2,pour_speed,g,density,E,nu,source_start,source_end);
            Add_Gravity(g);
            break;}

        case 42:{ // sand and thrown box
            particles.Store_Fp(true);
            Set_Grid(RANGE<TV>(TV(-.1,-.1,-.1),TV(3.1,1.1,1.1))*m,TV_INT(3,1,1));

            RIGID_BODY<TV>& bw=tests.Add_Analytic_Box(TV(3.2,0.1,1.2));
            bw.is_static=true;
            bw.Frame().t=TV(1.5,0.05,0.5);
            T h=0.3;
            RIGID_BODY<TV>& w1=tests.Add_Analytic_Box(TV(0.1,h,1.2));
            w1.Frame().t=TV(-0.05,h/2,0.5);
            w1.is_static=true;
            RIGID_BODY<TV>& w2=tests.Add_Analytic_Box(TV(0.1,h,1.2));
            w2.Frame().t=TV(3.05,h/2,0.5);
            w2.is_static=true;
            RIGID_BODY<TV>& w3=tests.Add_Analytic_Box(TV(3+0.2,h,0.1));
            w3.Frame().t=TV(1.5,h/2,-0.05);
            w3.is_static=true;
            RIGID_BODY<TV>& w4=tests.Add_Analytic_Box(TV(3+0.2,h,0.1));
            w4.Frame().t=TV(1.5,h/2,1.05);
            w4.is_static=true;

            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            RANGE<TV> sandbox(TV(0,0.1,0),TV(3,0.2,1));
            Seed_Particles(sandbox,0,0,density,particles_per_cell);
            LOG::printf("Particles: %d\n",particles.number);
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Add_Drucker_Prager(E,nu,foo_T2);
            Set_Lame_On_Particles(E,nu);
            TV g=m/(s*s)*TV(0,-9.81,0);
            Add_Gravity(g);
            RIGID_BODY<TV>& cube=tests.Add_Analytic_Box(TV(0.2,0.2,0.2),TV_INT()+1,density*2);
            cube.Frame().t=TV(0.3,0.6,0.5);
            cube.Frame().r=ROTATION<TV>::From_Euler_Angles((T)0.7,(T)0.7,(T)0.7);
            cube.Twist().linear=TV(4,-0.5,0);
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
            ARRAY_VIEW<VECTOR<T,3> >* colors=particles.template Get_Array<VECTOR<T,3> >("color");
            for(int p=0;p<particles.number;p++) (*colors)(p)=Sand_Color();
        } break;

        case 43:{ // Rigid sphere and sand pile
            particles.Store_Fp(true);
            Set_Grid(RANGE<TV>(TV(-.1,-.1,-.1),TV(1.1,1.1,1.1))*m);

            RIGID_BODY<TV>& bw=tests.Add_Analytic_Box(TV(1.2,0.1,1.2));
            bw.is_static=true;
            bw.Frame().t=TV(0.5,0.05,0.5);
            T h=0.2;
            RIGID_BODY<TV>& w1=tests.Add_Analytic_Box(TV(0.1,h,1.2));
            w1.Frame().t=TV(-0.05,h/2,0.5);
            w1.is_static=true;
            RIGID_BODY<TV>& w2=tests.Add_Analytic_Box(TV(0.1,h,1.2));
            w2.Frame().t=TV(1.05,h/2,0.5);
            w2.is_static=true;
            RIGID_BODY<TV>& w3=tests.Add_Analytic_Box(TV(1+0.2,h,0.1));
            w3.Frame().t=TV(0.5,h/2,-0.05);
            w3.is_static=true;
            RIGID_BODY<TV>& w4=tests.Add_Analytic_Box(TV(1+0.2,h,0.1));
            w4.Frame().t=TV(0.5,h/2,1.05);
            w4.is_static=true;

            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e6*unit_p*scale_E,nu=.3;
            RANGE<TV> sand(TV(0.3,0.1,0.3),TV(0.7,0.7,0.7));
            Seed_Particles(sand,0,0,density*1,particles_per_cell);
            LOG::printf("Particles: %d\n",particles.number);
            ARRAY_VIEW<VECTOR<T,3> >* colors=particles.template Get_Array<VECTOR<T,3> >("color");
            for(int p=0;p<particles.number;p++) (*colors)(p)=Sand_Color();
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Add_Drucker_Prager(E,nu,foo_T2);
            Set_Lame_On_Particles(E,nu);
            TV g=m/(s*s)*TV(0,-9.81,0);
            Add_Gravity(g);
            RIGID_BODY<TV>& sphere=tests.Add_Analytic_Sphere(0.075,density*10);
            sphere.Frame().t=TV(0.45,0.9,0.65);
            int sphere_id=sphere.particle_index;
            begin_time_step=[this,sphere_id](T time)
            {
                if(time<2)
                {
                    solid_body_collection.rigid_body_collection.rigid_body_particles.frame(sphere_id).t=TV(0.45,0.9,0.65);
                    solid_body_collection.rigid_body_collection.rigid_body_particles.twist(sphere_id).linear=TV(.1,0,0);
                }
            };
                
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
        } break;

        // ./mpm_rb -3d 160 -last_frame 20 -fooT3 1 -resolution 10 -fooT2 35 -use_exp_F -max_dt 2e-4 -symplectic_euler 
        case 160:{ // sand falling into a pile.
            particles.Store_Fp(true);
            Set_Grid(RANGE<TV>(TV(-.1,-0.02,-0.1),TV(0.1,0.06,0.1))*m,TV_INT(5,2,5));
            RANGE<TV> ground(TV(-10,-10,-10)*m,TV(10,0,10)*m);
            if(use_penalty_collisions) Add_Penalty_Collision_Object(ground);
            else Add_Collision_Object(ground,COLLISION_TYPE::stick,0);

            T E=35.37e4*unit_p*scale_E,nu=.3;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Add_Drucker_Prager(E,nu,foo_T2);
            TV spout(0,0.06*m,0);
            T density=(T)2200*unit_rho*scale_mass;
            T spout_width=8.334e-3*m;
            T pour_speed=.16319*m/s;
            TV gravity=TV(0,-9.8*m/(s*s),0);
            Set_Lame_On_Particles(E,nu);
            Add_Source(spout,TV(0,-1,0),spout_width/2,pour_speed,gravity,density,E,nu,0.08,foo_T3);
            Add_Gravity(gravity);
        } break;

        // ./mpm_rb -3d 161 -last_frame 20 -fooT3 1 -resolution 10 -fooT2 35 -use_exp_F -max_dt 2e-4 -grad_ls -rd -rd_stiffness 1e3
        case 161:{ // sand falling into a pile.
            particles.Store_Fp(true);
            Set_Grid(RANGE<TV>(TV(-.1,-0.02,-0.1),TV(0.1,0.06,0.1))*m,TV_INT(5,2,5));
            RIGID_BODY<TV>& ground=tests.Add_Analytic_Box(TV(10,0.1,10));
            ground.Frame().t=TV(0,-0.05,0);
            ground.is_static=true;
            T density=(T)2200*unit_rho*scale_mass;
            T E=35.37e4*unit_p*scale_E,nu=.3;
            T spout_width=8.334e-3*m;
            T pour_speed=.16319*m/s;
            TV gravity=TV(0,-9.8*m/(s*s),0);
            use_colored_sand=true;
            if(!no_implicit_plasticity) use_implicit_plasticity=true;
            Add_Drucker_Prager(E,nu,foo_T2);
            TV spout(0,0.06*m,0);
            Set_Lame_On_Particles(E,nu);
            Add_Source(spout,TV(0,-1,0),spout_width/2,pour_speed,gravity,density,E,nu,0.08,foo_T3);
            Add_Gravity(gravity);
        } break;
            
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
    if(forced_collision_type!=-1)
        for(int i=0;i<collision_objects.m;i++)
            collision_objects(i)->type=(COLLISION_TYPE)forced_collision_type;
}
template class STANDARD_TESTS<VECTOR<float,3> >;
template class STANDARD_TESTS<VECTOR<double,3> >;
}
