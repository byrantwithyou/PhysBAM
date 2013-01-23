//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORTEX_IMPLICIT_SURFACE
//#####################################################################
#ifndef __VORTEX_IMPLICIT_SURFACE__
#define __VORTEX_IMPLICIT_SURFACE__
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_IMPLICIT_SURFACE_LEVELSET.h>

namespace PhysBAM{
template<class TV>
struct VORTEX_IMPLICIT_SURFACE:public ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>
{
    typedef typename TV::SCALAR T;
    T k;
    virtual T f(const TV& X) const {return k-sin(X.x)*sin(X.y);}
    virtual TV df(const TV& X) const {return -TV(cos(X.x)*sin(X.y),sin(X.x)*cos(X.y));}
    virtual MATRIX<T,TV::m> ddf(const TV& X) const {T A=sin(X.x)*sin(X.y),B=-cos(X.x)*cos(X.y);return MATRIX<T,TV::m>(A,B,B,A);}
    virtual TV Closest_Point_Estimate(const TV& X) const {return (X-pi/2).Normalized()+pi/2;}
};
}
#endif
