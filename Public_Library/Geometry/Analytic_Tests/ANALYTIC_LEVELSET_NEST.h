//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_NEST
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_NEST__
#define __ANALYTIC_LEVELSET_NEST__
#include <Core/Arrays/ARRAY.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_NEST:public ANALYTIC_LEVELSET<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_LEVELSET<TV>* al;
    ARRAY<ANALYTIC_LEVELSET<TV>*> sub_al;

    ANALYTIC_LEVELSET_NEST(ANALYTIC_LEVELSET<TV>* ls): al(ls) {}
    ANALYTIC_LEVELSET_NEST* Add(ANALYTIC_LEVELSET<TV>* ls){sub_al.Append(ls);return this;}
    virtual ~ANALYTIC_LEVELSET_NEST();
    virtual T phi(const TV& X,T t,int& c) const;
    virtual TV N(const TV& X,T t,int c) const;
    virtual T dist(const TV& X,T t,int c) const;
    T find(const TV& X,T t,int c,TV* N) const;
};
}
#endif
