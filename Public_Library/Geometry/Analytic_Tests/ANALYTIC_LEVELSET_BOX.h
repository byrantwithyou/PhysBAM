//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_BOX
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_BOX__
#define __ANALYTIC_LEVELSET_BOX__
#include <Core/Math_Tools/RANGE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SIGNED.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_BOX:public ANALYTIC_LEVELSET_SIGNED<TV>
{
    typedef typename TV::SCALAR T;
    TV corner1,corner2;
    RANGE<TV> range;
    ANALYTIC_LEVELSET_BOX(TV corner1_in,TV corner2_in,int c_i,int c_o): ANALYTIC_LEVELSET_SIGNED<TV>(c_i,c_o),corner1(corner1_in),corner2(corner2_in),range(corner1_in,corner2_in) {}
    virtual T phi2(const TV& X,T t) const {return range.Signed_Distance(X);}
    virtual TV N2(const TV& X,T t) const {return range.Normal(X);}
};
}
#endif
