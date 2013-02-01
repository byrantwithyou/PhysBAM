//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_SIGNED
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_SIGNED__
#define __ANALYTIC_LEVELSET_SIGNED__
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_SIGNED:public ANALYTIC_LEVELSET<TV>
{
    typedef typename TV::SCALAR T;
    int color_i,color_o;
    ANALYTIC_LEVELSET_SIGNED(int color_i,int color_o);
    virtual ~ANALYTIC_LEVELSET_SIGNED() {}
    virtual T phi(const TV& X,T t,int& c) const;
    virtual TV N(const TV& X,T t,int c) const;
    virtual T dist(const TV& X,T t,int c) const;
    virtual T phi2(const TV& X,T t) const=0;
    virtual TV N2(const TV& X,T t) const=0;
};
}
#endif
