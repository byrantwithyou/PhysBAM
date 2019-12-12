//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_STATE__
#define __FLUID_STATE__
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV>
struct FLUID_STATE
{
    typedef typename TV::SCALAR T;
    FLUID_STATE()=default;
    virtual ~FLUID_STATE()=default;
};
}
#endif
