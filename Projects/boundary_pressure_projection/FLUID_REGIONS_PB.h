//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_REGIONS_PB__
#define __FLUID_REGIONS_PB__
#include "FLUID_REGIONS.h"
namespace PhysBAM{

template<class TV>
class FLUID_REGIONS_PB:public FLUID_REGIONS<TV>
{
public:
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    ARRAY<int,TV_INT> regions;

    FLUID_REGIONS_PB()=default;
    virtual ~FLUID_REGIONS_PB()=default;
};
}
#endif
