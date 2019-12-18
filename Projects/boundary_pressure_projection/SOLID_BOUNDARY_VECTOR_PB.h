//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SOLID_BOUNDARY_VECTOR_PB__
#define __SOLID_BOUNDARY_VECTOR_PB__
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Vectors/VECTOR.h>
#include "SOLID_BOUNDARY_VECTOR.h"
namespace PhysBAM{

template<class TV>
struct SOLID_BOUNDARY_VECTOR_PB:public SOLID_BOUNDARY_VECTOR<TV>
{
    typedef typename TV::SCALAR T;

    HASHTABLE<int,TV> V;
    HASHTABLE<int,TWIST<TV> > twist;

    SOLID_BOUNDARY_VECTOR_PB()=default;
    virtual ~SOLID_BOUNDARY_VECTOR_PB()=default;
};
}
#endif
