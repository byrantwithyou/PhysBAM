//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_VORTEX
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_VORTEX__
#define __ANALYTIC_LEVELSET_VORTEX__
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SIGNED.h>
#include <PhysBAM_Geometry/Analytic_Tests/VORTEX_IMPLICIT_SURFACE.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_VORTEX:public ANALYTIC_LEVELSET_SIGNED<TV>
{
    typedef typename TV::SCALAR T;
    T k;
    VORTEX_IMPLICIT_SURFACE<TV> vis;

    ANALYTIC_LEVELSET_VORTEX(T kk,int c_i,int c_o): ANALYTIC_LEVELSET_SIGNED<TV>(c_i,c_o),k(kk) {vis.k=k;}

    virtual T phi2(const TV& X,T t) const {return vis.Phi(X);}
    virtual TV N2(const TV& X,T t) const {return vis.Normal(X)*sign(vis.f(X));}
};
}
#endif
