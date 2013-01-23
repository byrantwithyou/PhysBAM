//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_IMPLICIT_SURFACE_LEVELSET
//#####################################################################
#ifndef __ANALYTIC_IMPLICIT_SURFACE_LEVELSET__
#define __ANALYTIC_IMPLICIT_SURFACE_LEVELSET__
#include <PhysBAM_Tools/Matrices/MATRIX.h>

namespace PhysBAM{
template<class TV>
class ANALYTIC_IMPLICIT_SURFACE_LEVELSET
{
    typedef typename TV::SCALAR T;
public:

    T tolerance;
    ANALYTIC_IMPLICIT_SURFACE_LEVELSET(): tolerance((T)1e-25) {}

    virtual T f(const TV& X) const=0;
    virtual TV df(const TV& X) const=0;
    virtual MATRIX<T,TV::m> ddf(const TV& X) const=0;
    virtual TV Closest_Point_Estimate(const TV& X) const {return X;}

    VECTOR<T,TV::m+1> Find_Closest_Point(const TV& X) const;
    T Phi(const TV& X) const;
    TV Normal(const TV& X) const;
};
}
#endif
