//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SOLID_BOUNDARY_VECTOR__
#define __SOLID_BOUNDARY_VECTOR__
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV>
struct SOLID_BOUNDARY_VECTOR
{
    typedef typename TV::SCALAR T;

    SOLID_BOUNDARY_VECTOR()=default;
    virtual ~SOLID_BOUNDARY_VECTOR()=default;
};
}
#endif
