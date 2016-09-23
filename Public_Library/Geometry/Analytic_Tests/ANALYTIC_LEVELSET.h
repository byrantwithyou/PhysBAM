//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET
//#####################################################################
#ifndef __ANALYTIC_LEVELSET__
#define __ANALYTIC_LEVELSET__
#include <Core/Vectors/VECTOR.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET
{
    typedef typename TV::SCALAR T;

    virtual ~ANALYTIC_LEVELSET();

    static T Large_Phi() {return 1000;}
    virtual T phi(const TV& X,T t,int& c) const=0;
    virtual TV N(const TV& X,T t,int c) const=0;
    virtual T dist(const TV& X,T t,int c) const=0; // signed distance to color
    void Test(const RANGE<TV>& domain) const;
};
}
#endif
