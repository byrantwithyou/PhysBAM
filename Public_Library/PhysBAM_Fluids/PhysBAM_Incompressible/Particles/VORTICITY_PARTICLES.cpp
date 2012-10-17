//#####################################################################
// Copyright 2004-2009, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORTICITY_PARTICLES
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> VORTICITY_PARTICLES<TV>::
VORTICITY_PARTICLES()
    :vorticity(0,0),radius(0,0)
{
    this->Add_Array(ATTRIBUTE_ID_VORTICITY,&vorticity);
    this->Add_Array(ATTRIBUTE_ID_RADIUS,&radius);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> VORTICITY_PARTICLES<TV>::
~VORTICITY_PARTICLES()
{}
//#####################################################################
template class VORTICITY_PARTICLES<VECTOR<float,1> >;
template class VORTICITY_PARTICLES<VECTOR<float,2> >;
template class VORTICITY_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VORTICITY_PARTICLES<VECTOR<double,1> >;
template class VORTICITY_PARTICLES<VECTOR<double,2> >;
template class VORTICITY_PARTICLES<VECTOR<double,3> >;
#endif
}
