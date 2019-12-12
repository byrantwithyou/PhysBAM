//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_BC__
#define __FLUID_BC__
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV>
class FLUID_BC
{
public:
    typedef typename TV::SCALAR T;
    FLUID_BC()=default;
    virtual ~FLUID_BC()=default;
};
}
#endif
