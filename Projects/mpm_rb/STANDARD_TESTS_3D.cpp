//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
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
#include <Deformables/Collisions_And_Interactions/PINNING_FORCE.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Constitutive_Models/MOONEY_RIVLIN_CURVATURE.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/OPENSUBDIV_SURFACE_CURVATURE_FORCE.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE_3D.h>
#include <Rigids/Forces_And_Torques/RIGID_GRAVITY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_SPHERE.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_DRUCKER_PRAGER.h>
#include <Hybrid_Methods/Forces/MPM_OLDROYD_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_VISCOSITY.h>
#include <Hybrid_Methods/Forces/VOLUME_PRESERVING_OB_NEO_HOOKEAN.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <fstream>
#include "POUR_SOURCE.h"
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

            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(TV((T).4,(T).5,(T).4)*m,TV((T).6,(T).7,(T).6)*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            for(int i=0;i<particles.number;i++){
                particles.X(i)=R.Rotate(particles.X(i)-c)+c;
                particles.V(i)=-t*vel;}

            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            TV g=m/(s*s)*TV(0,-1.8,0);
            Add_Gravity(g);
            Add_Collision_Object(Make_IO(PLANE<T>(n,c)));

            // Dump solution to viewer
            write_output_files=[=](int frame)
            {
                T time=frame*frame_dt;
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

            // ./mpm_rb -3d 1 -double -rd_stiffness 1e-1 -max_dt .01 -T 0.6 -rd_friction 0.1 -last_frame 200
        case 1:{
            T angle=extra_T(0);
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RIGID_BODY<TV>& ground=tests.Add_Analytic_Box(TV(1,0.2,1));
            ground.is_static=true;
            RIGID_BODY<TV>& cy1=tests.Add_Rigid_Body("skinnycyllink",(T).1*m,(T)0);
            cy1.Frame().t=TV(0,0.11,0)*m;
            cy1.Frame().r=ROTATION<TV>(pi/2,TV(1,0,0));
            RIGID_BODY<TV>& cy2=tests.Add_Rigid_Body("skinnycyllink",(T).1*m,(T)0);
            cy2.Frame().r=ROTATION<TV>(pi/2,TV(1,0,0))*ROTATION<TV>(angle,TV(0,0,1));
            cy2.Frame().t=(TV(0,0.13,0)+cy2.Frame().r.Rotate(TV(0,0.05,0)))*m;
            TV g=m/(s*s)*TV(0,-1.8,0);
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
        } break;

            // SIGGRAPH test 7, Wedging test, Rigid-MPM
            // sticking: ./mpm_rb -3d 71 -double -scale_E .1 -rd_stiffness 1e2 -max_dt .01 -regular_seeding -rd_friction 0.9
            // slipping: ./mpm_rb -3d 71 -double -scale_E .1 -rd_stiffness 1e2 -max_dt .01 -regular_seeding -rd_friction 0.3
        case 71:{
            Set_Grid(RANGE<TV>::Centered_Box()*2*m);
            RIGID_BODY<TV>& lw=tests.Add_Analytic_Box(TV(0.4,1,1));
            lw.Frame().t=TV(-0.2,0.5,0.5);
            lw.is_static=true;
            Add_Collision_Object(Make_IO(PLANE<T>(TV(-1,0,0),TV(1,0.5,0.5))));

            TV g=m/(s*s)*TV(0,-1.8,0);
            RIGID_BODY<TV>& cube=tests.Add_Analytic_Box(TV(0.4,0.4,0.4));
            cube.Frame().t=TV(0.2,0.5,0.5)*m;
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);

            SPHERE<TV> sphere(TV(.7,.5,0.5)*m,.4*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            Add_Gravity(g);
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
