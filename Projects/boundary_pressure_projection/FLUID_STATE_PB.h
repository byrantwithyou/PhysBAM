//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_STATE_PB__
#define __FLUID_STATE_PB__
#include "FLUID_STATE.h"
namespace PhysBAM{

template<class TV>
struct FLUID_STATE_PB:public FLUID_STATE<TV>
{
    typedef typename TV::SCALAR T;
    FLUID_STATE_PB()=default;
    virtual ~FLUID_STATE_PB()=default;
};
}
#endif
