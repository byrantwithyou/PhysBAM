//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_POISSON_SOLUTION
//#####################################################################
#ifndef __ANALYTIC_POISSON_SOLUTION__
#define __ANALYTIC_POISSON_SOLUTION__
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_POISSON_SOLUTION
{
    typedef typename TV::SCALAR T;

    virtual ~ANALYTIC_POISSON_SOLUTION(){}

    virtual T u(const TV& X)=0;
    virtual TV du(const TV& X)=0;
    virtual T f_volume(const TV& X)=0;
};
}
#endif
