//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SOLID_STATE_PB__
#define __SOLID_STATE_PB__
#include <Core/Vectors/VECTOR.h>
#include "SOLID_STATE.h"
namespace PhysBAM{

template<class TV>
struct SOLID_STATE_PB:public SOLID_STATE<TV>
{
    typedef typename TV::SCALAR T;

    T time;
    ARRAY<TV> X,V;
    ARRAY<FRAME<TV> > frame;
    ARRAY<TWIST<TV> > twist;

    SOLID_STATE_PB()=default;
    virtual ~SOLID_STATE_PB()=default;
};
}
#endif
