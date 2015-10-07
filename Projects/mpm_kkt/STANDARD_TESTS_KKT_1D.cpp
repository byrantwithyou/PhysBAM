//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Data_Structures/KD_TREE.h>
#include <Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Tools/Matrices/FRAME.h>
#include <Tools/Matrices/MATRIX.h>
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
            grid.Initialize(TV_INT()+resolution*2-1,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(TV(.45),TV(.55));
            T density=2*scale_mass;
            Seed_Particles_Helper(box,[=](const TV& X){return TV();},
                [=](const TV&){return MATRIX<T,1>();},
                density,particles_per_cell);
            T total_mass=particles.mass.Sum();
            TV total_momentum=particles.V.Weighted_Sum(particles.mass);
            TV dV=total_momentum/total_mass;
            particles.V-=dV;
            Add_Gravity(TV(-0.8));
        } break;
        case 2:{ // full box
            grid.Initialize(TV_INT()+resolution*2-1,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(grid.DX(),TV::All_Ones_Vector()-grid.DX());
            T density=scale_mass;
            Seed_Particles_Helper(box,0,0,density,particles_per_cell);
            Add_Gravity(TV(-1.8));
        } break; 
        case 3:{ // constant velocity
            grid.Initialize(TV_INT()+resolution*2-1,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(TV(.45),TV(.55));
            T density=2*scale_mass;
            Seed_Particles_Helper(box,[=](const TV& X){return TV(0.2);},
                [=](const TV&){return MATRIX<T,1>();},
                density,particles_per_cell);
        } break;
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
    // initialize coarse grid
    coarse_grid.Initialize(TV_INT()+resolution,grid.domain.Thickened(grid.dX/2),true);
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