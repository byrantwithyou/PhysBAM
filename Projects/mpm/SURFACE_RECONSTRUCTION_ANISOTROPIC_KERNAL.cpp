//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<TV>::
SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<TV>::
~SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL()
{}
//#####################################################################
template class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<VECTOR<float,2> >;
template class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<VECTOR<float,3> >;
template class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<VECTOR<double,2> >;
template class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<VECTOR<double,3> >;
}
