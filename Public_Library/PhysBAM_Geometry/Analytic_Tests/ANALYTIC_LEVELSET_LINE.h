//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_SIGNED
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_LINE__
#define __ANALYTIC_LEVELSET_LINE__
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SIGNED.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_LINE:public ANALYTIC_LEVELSET_SIGNED<TV>
{
    typedef typename TV::SCALAR T;
    typename BASIC_GEOMETRY_POLICY<TV>::HYPERPLANE plane; // N points outward
    ANALYTIC_LEVELSET_LINE(const TV& X,const TV& N,int c_i=0,int c_o=-4): ANALYTIC_LEVELSET_SIGNED<TV>(c_i,c_o),plane(N.Normalized(),X) {}
    virtual T phi2(const TV& X,T t) const {return plane.Signed_Distance(X);}
    virtual TV N2(const TV& X,T t) const {return plane.normal;}
};
}
#endif
