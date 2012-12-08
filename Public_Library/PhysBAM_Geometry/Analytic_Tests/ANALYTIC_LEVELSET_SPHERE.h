//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_SPHERE
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_SPHERE__
#define __ANALYTIC_LEVELSET_SPHERE__
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SIGNED.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_SPHERE:public ANALYTIC_LEVELSET_SIGNED<TV>
{
    typedef typename TV::SCALAR T;
    TV cen;
    T r;
    ANALYTIC_LEVELSET_SPHERE(TV cc,T rr,int c_i=0,int c_o=-4): ANALYTIC_LEVELSET_SIGNED<TV>(c_i,c_o),cen(cc),r(rr) {}
    virtual T phi2(const TV& X,T t) const {return (X-cen).Magnitude()-r;}
    virtual TV N2(const TV& X,T t) const {return (X-cen).Normalized();}
};
}
#endif
