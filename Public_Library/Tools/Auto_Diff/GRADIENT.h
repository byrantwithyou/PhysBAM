//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRADIENT
//##################################################################### 
#ifndef __GRADIENT__
#define __GRADIENT__

#include <Tools/Auto_Diff/VEC_HOLDER.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

template<class TV,class VEC_IN>
struct GRADIENT
{
    typedef typename TV::SCALAR T;
    typedef VEC_IN VEC;
    VEC x;

    GRADIENT operator+ () const {return *this;}
    GRADIENT<TV,decltype(VEC_NEG::Type(VEC()))> operator- () const {GRADIENT<TV,decltype(VEC_NEG::Type(VEC()))> r;VEC_NEG()(r.x,x);return r;}
    GRADIENT<TV,decltype(VEC_SCALE::Type(VEC(),T()))> operator* (T a) const {GRADIENT<TV,decltype(VEC_SCALE::Type(VEC(),T()))> r;VEC_SCALE()(r.x,x,a);return r;}
    GRADIENT<TV,decltype(VEC_SCALE_DIV::Type(VEC(),T()))> operator/ (T a) const {GRADIENT<TV,decltype(VEC_SCALE_DIV::Type(VEC(),T()))> r;VEC_SCALE_DIV()(r.x,x,a);return r;}

    template<class VEC1>
    GRADIENT<TV,decltype(VEC_ADD::Type(VEC(),VEC1()))> operator+ (const GRADIENT<TV,VEC1>& z) const {GRADIENT<TV,decltype(VEC_ADD::Type(VEC(),VEC1()))> r;VEC_ADD()(r.x,x,z.x);return r;}

    template<class VEC1>
    GRADIENT<TV,decltype(VEC_SUB::Type(VEC(),VEC1()))> operator- (const GRADIENT<TV,VEC1>& z) const {GRADIENT<TV,decltype(VEC_SUB::Type(VEC(),VEC1()))> r;VEC_SUB()(r.x,x,z.x);return r;}

    TV operator()(int i) const {TV v;Get(v,x,i);return v;}

    template<int i>
    void Set_Entry(const TV& v) {Set<i>(x,v);}
};

template<class TV,class VEC,class VEC1>
void Fill_From(GRADIENT<TV,VEC>& out,const GRADIENT<TV,VEC1>& in)
{Fill_From(out.x,in.x);}

template<class TV,class VEC>
GRADIENT<TV,decltype(VEC_SCALE::Type(VEC(),typename TV::SCALAR()))> operator* (typename TV::SCALAR a,const GRADIENT<TV,VEC>& u){return u*a;}

template<class TV,int d,class VEC> inline void
Extract(VECTOR<TV,d>& dx,const GRADIENT<TV,VEC>& v)
{Extract(dx,v.x);}

}
}
#endif
