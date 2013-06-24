//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header PARTICLES_FORWARD
//#####################################################################
#ifndef __PARTICLES_FORWARD__
#define __PARTICLES_FORWARD__
#include <Tools/Arrays/ATTRIBUTE_ID.h>

namespace PhysBAM{

template<class TV> class PARTICLES;

template<class T_PARTICLES> class PARTICLES_POOL;

template<class TV> class PARTICLES_CONNECTIVITY;
template<class TV,class T_PARTICLES> class PARTICLES_SUBSET;
const ATTRIBUTE_ID ATTRIBUTE_ID_X(1);
const ATTRIBUTE_ID ATTRIBUTE_ID_ID(20);
}
#endif
