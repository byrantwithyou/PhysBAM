//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/LAGGED_FORCE.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE_RB.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
#include <Hybrid_Methods/Forces/MPM_PLASTICITY_MODEL.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS_SPLINE.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_EXAMPLE_RB<TV>::
MPM_EXAMPLE_RB(const STREAM_TYPE stream_type)
    :MPM_EXAMPLE<TV>(stream_type)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_EXAMPLE_RB<TV>::
~MPM_EXAMPLE_RB()
{
}
//#####################################################################
namespace PhysBAM{
template class MPM_EXAMPLE_RB<VECTOR<float,2> >;
template class MPM_EXAMPLE_RB<VECTOR<float,3> >;
template class MPM_EXAMPLE_RB<VECTOR<double,2> >;
template class MPM_EXAMPLE_RB<VECTOR<double,3> >;
}
