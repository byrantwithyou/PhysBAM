//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SOLID_STATE__
#define __SOLID_STATE__
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV>
struct SOLID_STATE
{
    typedef typename TV::SCALAR T;

    SOLID_STATE()=default;
    virtual ~SOLID_STATE()=default;
};
}
#endif
