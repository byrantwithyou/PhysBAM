//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_LEVELSET_ROTATE
//#####################################################################
#ifndef __ANALYTIC_LEVELSET_ROTATE__
#define __ANALYTIC_LEVELSET_ROTATE__
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET.h>

namespace PhysBAM{
template<class TV>
struct ANALYTIC_LEVELSET_ROTATE:public ANALYTIC_LEVELSET<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_LEVELSET<TV>* al;
    T w;

    ANALYTIC_LEVELSET_ROTATE(ANALYTIC_LEVELSET<TV>* al,T w): al(al),w(w) {}
    ~ANALYTIC_LEVELSET_ROTATE() {delete al;}
    virtual T phi(const TV& X,T t,int& c) const {MATRIX<T,2> Q(MATRIX<T,2>::Rotation_Matrix(w*t));return al->phi(Q.Transpose_Times(X),t,c);}
    virtual TV N(const TV& X,T t,int c) const {MATRIX<T,2> Q(MATRIX<T,2>::Rotation_Matrix(w*t));return Q*(al->N(Q.Transpose_Times(X),t,c));}
    virtual T dist(const TV& X,T t,int c) const {MATRIX<T,2> Q(MATRIX<T,2>::Rotation_Matrix(w*t));return al->dist(Q.Transpose_Times(X),t,c);}
};
}
#endif
