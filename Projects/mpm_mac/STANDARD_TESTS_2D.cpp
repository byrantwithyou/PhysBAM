//#####################################################################
// Copyright 2015, Lin Huang, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/KD_TREE.h>
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
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Seeding/POISSON_DISK.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/PINNING_FORCE.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_GRAVITY.h>
#include <Hybrid_Methods/Forces/MPM_VISCOSITY.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <fstream>
#include "STANDARD_TESTS_2D.h"
namespace PhysBAM{
//#####################################################################
// Function Initialize_Implicit_Surface
//
// This was copied from DEFORMABLES_STANDARD_TESTS.cpp
// TODO: put this function somewhere more convenient (maybe as a constructor of LEVELSET_IMPLICIT_OBJECT?)
//#####################################################################
template<class T> LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >*
Initialize_Implicit_Surface(TRIANGULATED_SURFACE<T>& surface,int max_res)
{
    typedef VECTOR<int,3> TV_INT;
    LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >& undeformed_levelset=*LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >::Create();
    surface.Update_Bounding_Box();
    RANGE<VECTOR<T,3> > box=*surface.bounding_box;
    GRID<VECTOR<T,3> >& ls_grid=undeformed_levelset.levelset.grid;
    ARRAY<T,TV_INT>& phi=undeformed_levelset.levelset.phi;
    ls_grid=GRID<VECTOR<T,3> >::Create_Grid_Given_Cell_Size(box,box.Edge_Lengths().Max()/max_res,false,5);
    phi.Resize(ls_grid.Domain_Indices(3));
    LEVELSET_MAKER_UNIFORM<VECTOR<T,3> >::Compute_Level_Set(surface,ls_grid,3,phi);
    undeformed_levelset.Update_Box();
    return &undeformed_levelset;
}
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
    PHASE_ID number_phases(1);
    switch(test_number)
    {
        case 1:{ // stationary circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,0,0,density,particles_per_cell);
        } break;
        case 2:{ // translating circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(sphere,[=](const TV& X){return TV(m/s,0);},0,density,particles_per_cell);
        } break;
        case 3:{ // rotating circle
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            SPHERE<TV> sphere(TV(.5,.5)*m,.3*m);
            VECTOR<T,1> angular_velocity(0.4/s);
            T density=2*unit_rho*scale_mass;
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
            Seed_Particles(sphere,0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
            Add_Walls(-1,COLLISION_TYPE::slip,.1*m);
        } break;
        case 5:{ // stationary pool
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(RANGE<TV>(TV(.1*m,.1*m),TV(.9*m,.5*m)),0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
            Add_Walls(-1,COLLISION_TYPE::slip,.1*m);
        } break;
        case 6:{ // a dropping chunk of fluid
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            Seed_Particles(RANGE<TV>(TV(.2*m,.2*m),TV(.5*m,.8*m)),0,0,density,particles_per_cell);
            gravity=TV(0,-1)*m/sqr(s);
            Add_Walls(-1,COLLISION_TYPE::slip,.1*m);
        } break;
        case 7:{ // stationary circles in two phases
            number_phases=PHASE_ID(2);
            particles.Store_Phase(true);
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            SPHERE<TV> sphere0(TV(.3,.3)*m,.1*m);
            Seed_Particles(sphere0,0,0,density,particles_per_cell);
            int n=particles.phase.m;
            SPHERE<TV> sphere1(TV(.6,.6)*m,.1*m);
            Seed_Particles(sphere1,0,0,density,particles_per_cell);
            particles.phase.Array_View(n,particles.phase.m-n).Fill(1);
        } break;
        case 8:{ // concave shape
            this->use_phi=true;
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            SPHERE<TV> sphere0(TV(.3,.3)*m,.1*m);
            SPHERE<TV> sphere1(TV(.44,.44)*m,.1*m);
            typedef ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> > TOBJ;
            TOBJ* obj0=new TOBJ(sphere0);
            TOBJ* obj1=new TOBJ(sphere1);
            IMPLICIT_OBJECT_UNION<TV> shape(obj0,obj1);
            Seed_Particles(shape,0,0,density,particles_per_cell);
        } break;
        case 9:
        case 10:{ // freefall circles with different phases
            if(test_number==9) number_phases=PHASE_ID(2);
            particles.Store_Phase(true);
            Set_Grid(RANGE<TV>::Unit_Box()*m);
            T density=2*unit_rho*scale_mass;
            gravity=TV(0,-1)*m/sqr(s);
            SPHERE<TV> sphere0(TV(.3,.5)*m,.1*m);
            SPHERE<TV> sphere1(TV(.7,.5)*m,.1*m);
            Seed_Particles(sphere0,0,0,density,particles_per_cell);
            int n=particles.phase.m;
            Seed_Particles(sphere1,0,0,density,particles_per_cell);
            if(test_number==9) particles.phase.Array_View(n,particles.phase.m-n).Fill(1);
            Add_Walls(-1,COLLISION_TYPE::slip,.1*m);
            // a wall in the middle preventing the two circles from touching
            RANGE<TV> wall=RANGE<TV>(TV(.45,0)*m,TV(.55,1)*m);
            Add_Collision_Object(wall,COLLISION_TYPE::slip,0);
        } break;
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
    phases.Resize(number_phases);
    if(forced_collision_type!=-1)
        for(int i=0;i<collision_objects.m;i++)
            collision_objects(i)->type=(COLLISION_TYPE)forced_collision_type;
}
template class STANDARD_TESTS<VECTOR<float,2> >;
template class STANDARD_TESTS<VECTOR<double,2> >;
}
