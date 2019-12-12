//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_BC_PB__
#define __FLUID_BC_PB__
#include "FLUID_BC.h"
namespace PhysBAM{

template<class TV>
class FLUID_BC_PB:public FLUID_BC<TV>
{
public:
    typedef typename TV::SCALAR T;
    FLUID_BC_PB()=default;
    virtual ~FLUID_BC_PB()=default;
};
}
#endif
