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
#include "STANDARD_TESTS_2D.h"
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
    switch(test_number)
    {
        case 1:{ // stationary circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=unit_rho*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
        } break;
        case 2:{ // translating circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=unit_rho*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV(m/s,0);},0,density,particles_per_cell);
        } break;
        case 15:{ // colliding circles
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere0(TV(.25,.5)*m,.2*m);
            SPHERE<TV> sphere1(TV(.75,.5)*m,.2*m);
            T density=unit_rho*scale_mass;
            Seed_Particles(sphere0,[=](const TV& X){return TV(m/s,0);},0,density,particles_per_cell);
            Seed_Particles(sphere1,[=](const TV& X){return TV(-m/s,0);},0,density,particles_per_cell);
        } break;
        case 3:{ // rotating circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            VECTOR<T,1> angular_velocity(0.4/s);
            T density=unit_rho*scale_mass;
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
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            Add_Gravity(TV(0,-1)*m/sqr(s));
        } break;
        case 5:{ // stationary pool
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=unit_rho*scale_mass;
            Seed_Particles(RANGE<TV>(TV(),TV(m,.5*m)),0,0,density,particles_per_cell);
            Add_Gravity(TV(0,-1)*m/sqr(s));
        } break;
        case 6:{ // freefall rectangle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=unit_rho*scale_mass;
            Seed_Particles(RANGE<TV>(TV(.2*m,.2*m),TV(.5*m,.8*m)),0,0,density,particles_per_cell);
            Add_Gravity(TV(0,-1)*m/sqr(s));
        } break;
        case 7:{ // (fluid test) dam break
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(grid.dX*2,TV(0.2,0.75)*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Quasi_Pressure(1e3*unit_p*scale_E,7);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 8:{ // (fluid test) circle drop
            // one: ./mpm -kkt -scale_E 0
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.7)*m,.2*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            Add_Quasi_Pressure(1e3*unit_p*scale_E,7);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 9:{ // (fluid test) pool of water w/ single particle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(grid.dX*(T).5,TV(1*m-grid.dX(0)*(T).5,0.25*m));
            T density=2*unit_rho*scale_mass;
            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Particle(TV(.5,.9),0,0,mass,volume);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 10:{ // Rayleigh Taylor
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(grid.dX*2,TV(1*m-2*grid.dX(0),0.20*m));
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            density*=10;
            box+=TV(0,0.20*m-2*grid.dX(0));
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 11:{ // full box
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            RANGE<TV> box(grid.dX*(T)2,TV::All_Ones_Vector()*m-grid.dX*2);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 12:{ // single particle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            T volume=grid.dX.Product()/particles_per_cell;
            T mass=density*volume;
            Add_Particle(TV(.8,.5),0,0,mass,volume);
            Add_Gravity(m/(s*s)*TV(0,-1.8));
        } break;
        case 13:{ // (fluid test) dam break; Rabecca Brannon test
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
            Add_Gravity(m/(s*s)*TV(0,-9.81));
            end_time_step=[=](T time)
                {
                    for(int i=0;i<particles.F.m;i++){
                        T J=particles.F(i).Determinant();
                        particles.F(i)=MATRIX<T,TV::m>()+pow<1,TV::m>(J);}
                    if(dt<1e-7) PHYSBAM_FATAL_ERROR("DT TOO SMALL");
                    for(int i=0;i<particles.V.m;i++)
                        if(particles.V(i).Magnitude()>50)
                            PHYSBAM_FATAL_ERROR("VELOCITY TOO BIG");
                };
        } break;

        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
}
template class STANDARD_TESTS<VECTOR<float,2> >;
template class STANDARD_TESTS<VECTOR<double,2> >;
}
