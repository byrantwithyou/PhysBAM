//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SOLID_BC_PB__
#define __SOLID_BC_PB__
#include <Core/Vectors/VECTOR.h>
#include "SOLID_BC.h"
namespace PhysBAM{

template<class TV>
class SOLID_BC_PB:public SOLID_BC<TV>
{
public:
    typedef typename TV::SCALAR T;

    SOLID_BC_PB()=default;
    virtual ~SOLID_BC_PB()=default;
};
}
#endif
