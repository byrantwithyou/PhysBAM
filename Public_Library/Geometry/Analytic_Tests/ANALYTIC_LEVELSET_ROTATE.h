//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_ROTATE
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_ROTATE__
#define __ANALYTIC_LEVELSET_ROTATE__
#include <Core/Matrices/ROTATION.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_ROTATE:public ANALYTIC_LEVELSET<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_LEVELSET<TV>* al;
    typename TV::SPIN w;
    TV center;

    ANALYTIC_LEVELSET_ROTATE(ANALYTIC_LEVELSET<TV>* al,const typename TV::SPIN& w,TV cc=TV()): al(al),w(w),center(cc) {}
    ~ANALYTIC_LEVELSET_ROTATE() {delete al;}
    virtual T phi(const TV& X,T t,int& c) const {ROTATION<TV> rot(ROTATION<TV>::From_Rotation_Vector(w*t));return al->phi(rot.Inverse_Rotate(X-center)+center,t,c);}
    virtual TV N(const TV& X,T t,int c) const {ROTATION<TV> rot(ROTATION<TV>::From_Rotation_Vector(w*t));return rot.Rotate(al->N(rot.Inverse_Rotate(X-center)+center,t,c));}
    virtual T dist(const TV& X,T t,int c) const {ROTATION<TV> rot(ROTATION<TV>::From_Rotation_Vector(w*t));return al->dist(rot.Inverse_Rotate(X-center)+center,t,c);}
};
}
#endif
