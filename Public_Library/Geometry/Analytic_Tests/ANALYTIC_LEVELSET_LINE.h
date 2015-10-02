//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_LINE
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_LINE__
#define __ANALYTIC_LEVELSET_LINE__
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SIGNED.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_LINE:public ANALYTIC_LEVELSET_SIGNED<TV>
{
    typedef typename TV::SCALAR T;
    TV x,n;
    ANALYTIC_LEVELSET_LINE(const TV& X,const TV& N,int c_i,int c_o): ANALYTIC_LEVELSET_SIGNED<TV>(c_i,c_o),x(X),n(N.Normalized()) {}
    virtual T phi2(const TV& X,T t) const {return n.Dot(X-x);}
    virtual TV N2(const TV& X,T t) const {return n;}
};
}
#endif
