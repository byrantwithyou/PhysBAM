//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_SCALE
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_SCALE__
#define __ANALYTIC_LEVELSET_SCALE__
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_SCALE:public ANALYTIC_LEVELSET<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_LEVELSET<TV>* al;
    T scale;

    ANALYTIC_LEVELSET_SCALE(ANALYTIC_LEVELSET<TV>* al,T scale): al(al),scale(scale) {}
    ~ANALYTIC_LEVELSET_SCALE() {delete al;}
    virtual T phi(const TV& X,T t,int& c) const {return (1+t*scale)*al->phi(X/(1+t*scale),t,c);}
    virtual TV N(const TV& X,T t,int c) const {return al->N(X/(1+t*scale),t,c);}
    virtual T dist(const TV& X,T t,int c) const {return (1+t*scale)*al->dist(X/(1+t*scale),t,c);}
};
}
#endif





