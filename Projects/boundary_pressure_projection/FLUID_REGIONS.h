//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_REGIONS__
#define __FLUID_REGIONS__
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV>
class FLUID_REGIONS
{
public:
    typedef typename TV::SCALAR T;

    int num_regions=0;

    FLUID_REGIONS()=default;
    virtual ~FLUID_REGIONS()=default;
};
}
#endif
