//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_TRANSLATE
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_TRANSLATE__
#define __ANALYTIC_LEVELSET_TRANSLATE__
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_TRANSLATE:public ANALYTIC_LEVELSET<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_LEVELSET<TV>* al;
    TV vel;

    ANALYTIC_LEVELSET_TRANSLATE(ANALYTIC_LEVELSET<TV>* al,const TV& vel): al(al),vel(vel) {}
    ~ANALYTIC_LEVELSET_TRANSLATE() {delete al;}
    virtual T phi(const TV& X,T t,int& c) const {return al->phi(X-vel*t,t,c);}
    virtual TV N(const TV& X,T t,int c) const {return al->N(X-vel*t,t,c);}
    virtual T dist(const TV& X,T t,int c) const {return al->dist(X-vel*t,t,c);}
};
}
#endif
