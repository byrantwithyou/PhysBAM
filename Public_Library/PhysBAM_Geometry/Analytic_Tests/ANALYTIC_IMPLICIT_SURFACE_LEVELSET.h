//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_IMPLICIT_SURFACE_LEVELSET
//#####################################################################
#ifndef __ANALYTIC_IMPLICIT_SURFACE_LEVELSET__
#define __ANALYTIC_IMPLICIT_SURFACE_LEVELSET__
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SIGNED.h>

namespace PhysBAM{
template<class TV>
class ANALYTIC_IMPLICIT_SURFACE_LEVELSET:public ANALYTIC_LEVELSET_SIGNED<TV>
{
    typedef typename TV::SCALAR T;
public:

    T tolerance;
    ANALYTIC_IMPLICIT_SURFACE_LEVELSET(int color_i,int color_o): ANALYTIC_LEVELSET_SIGNED<TV>(color_i,color_o),tolerance((T)1e-25) {}
    virtual ~ANALYTIC_IMPLICIT_SURFACE_LEVELSET(){}

    virtual T f(const TV& X,T t) const=0;
    virtual TV df(const TV& X,T t) const=0;
    virtual MATRIX<T,TV::m> ddf(const TV& X,T t) const=0;
    virtual TV Closest_Point_Estimate(const TV& X,T t) const {return X;}

    VECTOR<T,TV::m+1> Find_Closest_Point(const TV& X,T t) const;
    virtual T phi2(const TV& X,T t) const;
    virtual TV N2(const TV& X,T t) const;
};
}
#endif
