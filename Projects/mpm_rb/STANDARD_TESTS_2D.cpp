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
#include <Geometry/Basic_Geometry/LINE_2D.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/PLANE.h>
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
#include <Rigids/Forces_And_Torques/RIGID_GRAVITY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/PINNING_FORCE.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Constitutive_Models/QUASI_INCOMPRESSIBLE_FORCE.h>
#include <Deformables/Constitutive_Models/TAIT_PRESSURE_FORCE.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_OBJECT.h>
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
#include <fstream>
#include "POUR_SOURCE.h"
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
// Function Write_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,2> >::
Write_Output_Files(const int frame)
{
    if(write_output_files) write_output_files(frame);
    BASE::Write_Output_Files(frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,2> >::
Read_Output_Files(const int frame)
{
    if(read_output_files) read_output_files(frame);
    BASE::Read_Output_Files(frame);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void STANDARD_TESTS<VECTOR<T,2> >::
Initialize()
{
    this->asymmetric_system=true;
    switch(test_number)
    {
        case 1:{ // half-full box
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(TV(),TV(1,(T).5)*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;

            // ./mpm_rb 2 -float -sound_cfl -strong_cfl -reflection_bc -1 -symplectic_euler -scale_E .1 -coll_pair
        case 2:{ // half-full box with rigid body
            PHYSBAM_ASSERT(sizeof(T)==sizeof(float));
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(TV(),TV(1,(T).5)*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            TV g=m/(s*s)*TV(0,-1.8);
            Add_Gravity(g);
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("circle",(T).2*m,(T).5);
            rigid_body.Frame().t=TV((T)0.5,(T)0.75)*m;
            rigid_body.Set_Mass(2*kg);
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
        } break;

            // ./mpm_rb 3 -float -symplectic_euler -coll_pair
        case 3:{ // rigid sphere on rigid ground
            PHYSBAM_ASSERT(sizeof(T)==sizeof(float));
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            tests.Add_Ground((T).5,(T).1,0,1);
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("circle",(T).2*m,(T).5);
            rigid_body.Frame().t=TV((T)0.5,(T)0.75)*m;
            rigid_body.Set_Mass(1*kg);
            TV g=m/(s*s)*TV(0,-1.8);
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
        } break;

            // ./mpm_rb 4 -float -symplectic_euler -coll_pair
        case 4:{ // rigid sphere on rigid ground, bouncy
            PHYSBAM_ASSERT(sizeof(T)==sizeof(float));
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            tests.Add_Ground((T).5,(T).1,1,1);
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("circle",(T).2*m,(T).5);
            rigid_body.Frame().t=TV((T)0.5,(T)0.75)*m;
            rigid_body.Set_Mass(1*kg);
            rigid_body.coefficient_of_restitution=1;
            TV g=m/(s*s)*TV(0,-1.8);
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
        } break;

            // ./mpm_rb 3 -float -symplectic_euler -coll_pair
        case 5:{ // Many rigid spheres in box
            PHYSBAM_ASSERT(sizeof(T)==sizeof(float));
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RIGID_BODY<TV>& wl=tests.Add_Analytic_Box(TV(1,1));
            wl.Frame().t=TV(-.5,.5);
            wl.is_static=true;
            RIGID_BODY<TV>& wb=tests.Add_Analytic_Box(TV(1,1));
            wb.Frame().t=TV(.5,-.5);
            wb.is_static=true;
            RIGID_BODY<TV>& wr=tests.Add_Analytic_Box(TV(1,1));
            wr.Frame().t=TV(1.5,.5);
            wr.is_static=true;
            for(int i=0;i<10;i++){
                RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("circle",(T).2,(T).5);
                rigid_body.Frame().t=TV((T)0.5,(T)0.75+i*.5);}
            TV g=m/(s*s)*TV(0,-1.8);
            auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            solid_body_collection.rigid_body_collection.Add_Force(rg);
        } break;

            // ./mpm_rb 6 -float -scale_E .1  -last_frame 20 -rd_stiffness 1e3
        case 6:{ // MPM in box
            PHYSBAM_ASSERT(sizeof(T)==sizeof(float));
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(TV((T).2,(T).2)*m,TV((T).5,(T).5)*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            TV g=m/(s*s)*TV(0,-1.8);
            Add_Gravity(g);
            auto* wl=new ANALYTIC_IMPLICIT_OBJECT<LINE_2D<T> >(LINE_2D<T>(TV(1,0),TV(0.1,0.1)));
            auto* wb=new ANALYTIC_IMPLICIT_OBJECT<LINE_2D<T> >(LINE_2D<T>(TV(0,1),TV(0.1,0.1)));
            auto* wr=new ANALYTIC_IMPLICIT_OBJECT<LINE_2D<T> >(LINE_2D<T>(TV(-1,0),TV(.9,.1)));
            Add_Collision_Object(wl);
            Add_Collision_Object(wb);
            Add_Collision_Object(wr);
        } break;

            // ./mpm_rb 7 -float -scale_E .1 -rd_stiffness 1e3
        case 7:{ // MPM block and rigid circle, no collisions with boundary.
            PHYSBAM_ASSERT(sizeof(T)==sizeof(float));
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(TV(.2,.2)*m,TV(.8,(T).5)*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            // TV g=m/(s*s)*TV(0,-1.8);
            // Add_Gravity(g);
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("circle",(T).2*m,(T).5);
            rigid_body.Frame().t=TV((T)0.5,(T)0.75)*m;
            rigid_body.Set_Mass(2*kg);
            rigid_body.Twist().linear=TV(0,-1)*m/s;
            // auto* rg=new RIGID_GRAVITY<TV>(solid_body_collection.rigid_body_collection,0,g);
            // solid_body_collection.rigid_body_collection.Add_Force(rg);
        } break;

            // sticking ./mpm_rb 8 -float -scale_E .1 -rd_stiffness 1e1 -max_dt .01 -T .1 -T .1 -regular_seeding
            // sliding: ./mpm_rb 8 -float -scale_E .1 -rd_stiffness 1e1 -max_dt .01   -T .1 -T .1 -regular_seeding -rd_friction 0.1
        case 8:{ // MPM block and rigid circle, no collisions with boundary.
            T angle=extra_T(0);
            T vel=extra_T(1);
            PHYSBAM_ASSERT(angle>=0 && angle<=pi/2 && vel>=0);
            TV n(-sin(angle),cos(angle));
            TV t(cos(angle),sin(angle));
            MATRIX<T,2> M(t,n);
            TV c(.5,.5);
            
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(TV((T).4,(T).5)*m,TV((T).6,(T).7)*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            for(int i=0;i<particles.number;i++){
                particles.X(i)=M*(particles.X(i)-c)+c;
                particles.V(i)=-t*vel;}

            Add_Fixed_Corotated(1e3*unit_p*scale_E,0.3);
            TV g=m/(s*s)*TV(0,-1.8);
            Add_Gravity(g);
            Add_Collision_Object(Make_IO(LINE_2D<T>(n,c)));

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

        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
    if(forced_collision_type!=-1)
        for(int i=0;i<collision_objects.m;i++)
            collision_objects(i)->type=(COLLISION_TYPE)forced_collision_type;
}
template class STANDARD_TESTS<VECTOR<float,2> >;
template class STANDARD_TESTS<VECTOR<double,2> >;
}
