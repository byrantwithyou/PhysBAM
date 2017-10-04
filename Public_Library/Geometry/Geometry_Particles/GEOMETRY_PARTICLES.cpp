//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Vectors/VECTOR.h>
#include <Tools/Particles/PARTICLES.h>
#include <Tools/Particles/PARTICLES_FORWARD.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GEOMETRY_PARTICLES<TV>::
GEOMETRY_PARTICLES()
    :X(0,0),V(0,0),store_velocity(false)
{
    Add_Array("X",&X);
}
//#####################################################################
// Destructor 
//#####################################################################
template<class TV> GEOMETRY_PARTICLES<TV>::
~GEOMETRY_PARTICLES()
{}
//#####################################################################
template class GEOMETRY_PARTICLES<VECTOR<float,1> >;
template class GEOMETRY_PARTICLES<VECTOR<float,2> >;
template class GEOMETRY_PARTICLES<VECTOR<float,3> >;
template class GEOMETRY_PARTICLES<VECTOR<double,1> >;
template class GEOMETRY_PARTICLES<VECTOR<double,2> >;
template class GEOMETRY_PARTICLES<VECTOR<double,3> >;
}
