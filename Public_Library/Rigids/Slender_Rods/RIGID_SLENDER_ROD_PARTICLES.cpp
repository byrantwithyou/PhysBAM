//#####################################################################
// Copyright 2012, Steve Cook, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Rigids/Slender_Rods/RIGID_SLENDER_ROD_PARTICLES.h>

namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_SLENDER_ROD_PARTICLES<TV>::
RIGID_SLENDER_ROD_PARTICLES()
{
    Add_Array("mass",&mass);
    Add_Array("twist",&twist);
    Add_Array("X",&X);
    Add_Array("orientation",&orientation);
    Add_Array("length",&length);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_SLENDER_ROD_PARTICLES<TV>::
~RIGID_SLENDER_ROD_PARTICLES()
{}
template class RIGID_SLENDER_ROD_PARTICLES<VECTOR<float,1> >;
template class RIGID_SLENDER_ROD_PARTICLES<VECTOR<float,2> >;
template class RIGID_SLENDER_ROD_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_SLENDER_ROD_PARTICLES<VECTOR<double,1> >;
template class RIGID_SLENDER_ROD_PARTICLES<VECTOR<double,2> >;
template class RIGID_SLENDER_ROD_PARTICLES<VECTOR<double,3> >;
#endif
}
