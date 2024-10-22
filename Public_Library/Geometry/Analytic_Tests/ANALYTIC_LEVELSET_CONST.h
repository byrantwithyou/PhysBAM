//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_CONST
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_CONST__
#define __ANALYTIC_LEVELSET_CONST__
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SIGNED.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_CONST:public ANALYTIC_LEVELSET_SIGNED<TV>
{
    typedef typename TV::SCALAR T;
    T value;
    ANALYTIC_LEVELSET_CONST(T value,int c_i,int c_o)
        :ANALYTIC_LEVELSET_SIGNED<TV>(c_i,c_o),value(value)
    {}
    virtual T phi2(const TV& X,T t) const {return value;}
    virtual TV N2(const TV& X,T t) const {return TV::Axis_Vector(0);}
};
}
#endif
