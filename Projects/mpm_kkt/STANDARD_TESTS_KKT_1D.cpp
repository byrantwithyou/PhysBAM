//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/KD_TREE.h>
#include <Core/Math_Tools/RANGE_ITERATOR.h>
#include <Core/Matrices/FRAME.h>
#include <Core/Matrices/MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/PINNING_FORCE.h>
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
#include "STANDARD_TESTS_KKT_1D.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS_KKT<VECTOR<T,1> >::
STANDARD_TESTS_KKT(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_KKT_BASE<TV>(stream_type_input,parse_args),
    use_surface_tension(false),Nsurface(0)
{
    parse_args.Parse();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> STANDARD_TESTS_KKT<VECTOR<T,1> >::
~STANDARD_TESTS_KKT()
{
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,1> >::
Write_Output_Files(const int frame)
{
    BASE::Write_Output_Files(frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,1> >::
Read_Output_Files(const int frame)
{
    BASE::Read_Output_Files(frame);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,1> >::
Initialize()
{
    switch(test_number)
    {
        case 1:{ // circle free fall
            Set_Grid(RANGE<TV>::Unit_Box());
            RANGE<TV> box(TV(.45),TV(.55));
            T density=2*scale_mass;
            Seed_Particles(box,[=](const TV& X){return TV();},
                [=](const TV&){return MATRIX<T,1>();},
                density,particles_per_cell);
            T total_mass=particles.mass.Sum();
            TV total_momentum=particles.V.Weighted_Sum(particles.mass);
            TV dV=total_momentum/total_mass;
            particles.V-=dV;
            Add_Gravity(TV(-0.8));
            Add_Fluid_Wall(new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(RANGE<TV>(TV(-5),TV(0))));
            Add_Fluid_Wall(new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(RANGE<TV>(TV(1),TV(5))));
        } break;
        case 2:{ // full box
            Set_Grid(RANGE<TV>::Unit_Box());
            RANGE<TV> box(grid.dX,TV::All_Ones_Vector()-grid.dX);
            T density=scale_mass;
            Seed_Particles(box,0,0,density,particles_per_cell);
            Add_Gravity(TV(-1.8));
        } break; 
        case 3:{ // constant velocity
            Set_Grid(RANGE<TV>::Unit_Box());
            RANGE<TV> box(TV(.45),TV(.55));
            T density=2*scale_mass;
            Seed_Particles(box,[=](const TV& X){return TV(0.2);},
                [=](const TV&){return MATRIX<T,1>();},
                density,particles_per_cell);
        } break;
        case 4:{ // constant velocity
            Set_Grid(RANGE<TV>::Unit_Box());
            RANGE<TV> box(TV(.45),TV(.55));
            T density=scale_mass;
            Seed_Particles(box,[=](const TV& X){return TV(0.2);},
                [=](const TV&){return MATRIX<T,1>();},
                density,particles_per_cell);
            particles.lambda.Fill(1);
        } break;
        case 5:{ // gravity with one_over_lambda
            Set_Grid(RANGE<TV>::Unit_Box());
            RANGE<TV> box(TV(.45),TV(.55));
            T density=scale_mass;
            Seed_Particles(box,[=](const TV& X){return TV(0.2);},0,
                density,particles_per_cell);
            particles.lambda.Fill(1);
            Add_Gravity(TV(-0.8));
        } break;
        case 6:{ // analytic test
            Set_Grid(RANGE<TV>::Unit_Box());
            RANGE<TV> box(TV(.4),TV(.6));
            T density=scale_mass;
            Seed_Particles(box,[=](const TV& X){return TV(0.5*cos((T)pi/0.2*(X(0)-0.4)));},
                [=](const TV&){return MATRIX<T,1>();},
                density,particles_per_cell);
            particles.lambda.Fill(1);
        } break;
        case 7:{ // half box
            Set_Grid(RANGE<TV>::Unit_Box());
            RANGE<TV> box(TV(0.5),TV(1));
            T density=scale_mass;
            Seed_Particles(box,[=](const TV& X){return TV(0.5);},0,density,particles_per_cell);
            Add_Fluid_Wall(new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(RANGE<TV>(TV(-5),TV(0))));
            Add_Fluid_Wall(new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(RANGE<TV>(TV(1),TV(5))));
        } break;
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
}
//#####################################################################
// Function Begin_Frame
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,1> >::
Begin_Frame(const int frame)
{
}
//#####################################################################
// Function End_Frame
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,1> >::
End_Frame(const int frame)
{
}
//#####################################################################
// Function Begin_Time_Step
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,1> >::
Begin_Time_Step(const T time)
{
}
//#####################################################################
// Function End_Time_Step
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,1> >::
End_Time_Step(const T time)
{
}
template class STANDARD_TESTS_KKT<VECTOR<float,1> >;
template class STANDARD_TESTS_KKT<VECTOR<double,1> >;
}
