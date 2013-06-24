//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_ELLIPSOID
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_ELLIPSOID__
#define __ANALYTIC_LEVELSET_ELLIPSOID__
#include <Geometry/Analytic_Tests/ANALYTIC_IMPLICIT_SURFACE_LEVELSET.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_ELLIPSOID:public ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>
{
    typedef typename TV::SCALAR T;
    TV cen;
    TV r;
    MATRIX<T,TV::m> jacobian;
    ANALYTIC_LEVELSET_ELLIPSOID(TV cc,TV r_in,int c_i,int c_o): ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>(c_i,c_o),cen(cc),r(r_in) {
        for(int i=0;i<TV::m;i++) jacobian(i,i)=2*r.Magnitude_Squared()/r(i)/r(i);
    }
    virtual T f(const TV& X,T t) const {return (((X-cen)*(X-cen)/r/r).Sum()-(T)1)*r.Magnitude_Squared();}
    virtual TV df(const TV& X,T t) const {return (X-cen)/r/r*(T)2*r.Magnitude_Squared();}
    virtual MATRIX<T,TV::m> ddf(const TV& X,T t) const {return jacobian;}
    virtual TV Closest_Point_Estimate(const TV& X,T t) const {return X;}
};
}
#endif
