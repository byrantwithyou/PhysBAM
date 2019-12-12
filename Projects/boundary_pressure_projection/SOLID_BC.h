//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SOLID_BC__
#define __SOLID_BC__
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV>
class SOLID_BC
{
public:
    typedef typename TV::SCALAR T;

    SOLID_BC()=default;
    virtual ~SOLID_BC()=default;
};
}
#endif
