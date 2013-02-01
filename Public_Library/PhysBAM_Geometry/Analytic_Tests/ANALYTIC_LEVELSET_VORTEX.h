//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_VORTEX
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_VORTEX__
#define __ANALYTIC_LEVELSET_VORTEX__
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_IMPLICIT_SURFACE_LEVELSET.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_VORTEX:public ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>
{
    typedef typename TV::SCALAR T;
    T k;

    ANALYTIC_LEVELSET_VORTEX(T kk,int c_i,int c_o): ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>(c_i,c_o),k(kk) {}
    virtual ~ANALYTIC_LEVELSET_VORTEX() {}

    virtual T f(const TV& X,T t) const {return k-sin(X.x)*sin(X.y);}
    virtual TV df(const TV& X,T t) const {return -TV(cos(X.x)*sin(X.y),sin(X.x)*cos(X.y));}
    virtual MATRIX<T,TV::m> ddf(const TV& X,T t) const {T A=sin(X.x)*sin(X.y),B=-cos(X.x)*cos(X.y);return MATRIX<T,TV::m>(A,B,B,A);}
    virtual TV Closest_Point_Estimate(const TV& X,T t) const {return (X-pi/2).Normalized()+pi/2;}
};
}
#endif
