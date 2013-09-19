//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRADIENT_VEC
//##################################################################### 
#ifndef __GRADIENT_VEC__
#define __GRADIENT_VEC__

#include <Tools/Auto_Diff/COL_MAT.h>
#include <Tools/Auto_Diff/GRADIENT.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

#define MK_GV(TV,FLAGS) GRADIENT_VEC<TV,decltype(FLAGS)>

template<class TV,class H>
struct GRADIENT_VEC
{
    typedef typename TV::SCALAR T;
    // Note that this represents the result as a (n d) x d matrix (block column vector)
    // rather than as a block row vector, so the result is really transposed.
    COL_MAT<TV,H> c;

    MK_GV(TV,-H()) operator- () const {MK_GV(TV,-H()) r;r.c.Neg(c);return r;}
    GRADIENT_VEC operator+ () const {return *this;}
    MK_GV(TV,H()*1) operator* (T a) const {MK_GV(TV,H()*1) r;r.c.Scale(c,a);return r;}
    MK_GV(TV,H()/1) operator/ (T a) const {return *this*(1/a);}

    template<class H2>
    MK_GV(TV,H()+H2()) operator+ (const GRADIENT_VEC<TV,H2>& z) const
    {MK_GV(TV,H()+H2()) r;r.c.Add(c,z.c);return r;}

    template<class H2>
    MK_GV(TV,H()-H2()) operator- (const GRADIENT_VEC<TV,H2>& z) const
    {MK_GV(TV,H()-H2()) r;r.c.Sub(c,z.c);return r;}

    GRADIENT<TV,CV_FLAGS<H::n,0> > Transpose_Times(ZERO_VEC<TV> v) const
    {return GRADIENT<TV,CV_FLAGS<H::n,0> >();}

    GRADIENT<TV,CV_FLAGS<H::n,H::mask|H::sp_mask> > Transpose_Times(const TV& v) const
    {
        GRADIENT<TV,CV_FLAGS<H::n,H::mask|H::sp_mask> > g;
        Times_Helper(g.c,c,v);
        return g;
    }

    template<int i>
    void Set_Entry(const TV& v) {STATIC_ASSERT(H::mask&(1<<i));c.Set_Entry(v,(VECTOR<int,i>*)0);}

    template<class H2>
    void Fill_From(const GRADIENT_VEC<TV,H2>& z)
    {STATIC_ASSERT_SAME(H,decltype(Choose(H(),H2())));c.Fill_From(z.c);}
};

template<class TV,class H>
MK_GV(TV,H()*1) operator* (typename TV::SCALAR a,const GRADIENT_VEC<TV,H>& u)
{return u*a;}

template<int m,int n,int m2> CM_FLAGS<n,m2*m,0> Outer(V_FLAGS<m>,CV_FLAGS<n,m2>);

template<class TV,class H> GRADIENT_VEC<TV,CM_FLAGS<H::n,H::mask,0> >
Outer_Product(const TV& u,const GRADIENT<TV,H>& v)
{
    GRADIENT_VEC<TV,CM_FLAGS<H::n,H::mask,0> > GV;
    Outer_Product_Helper(GV.c,v.c,u);
    return GV;
}

template<class TV,class H> GRADIENT_VEC<TV,CM_FLAGS<H::n,0,0> >
Outer_Product(const ZERO_VEC<TV>& u,const GRADIENT<TV,H>& v)
{return GRADIENT_VEC<TV,CM_FLAGS<H::n,0,0> >();}

template<class TV,class H>
GRADIENT_VEC<TV,CM_FLAGS<H::n,0,0> > operator* (const ZERO_MAT<TV>& a,const GRADIENT_VEC<TV,H>& u){return GRADIENT_VEC<TV,CM_FLAGS<H::n,0,0> >();}

template<class TV,class H>
GRADIENT_VEC<TV,H> operator* (const ID_MAT<TV>& a,const GRADIENT_VEC<TV,H>& u){return u;}

template<class TV,class H>
MK_GV(TV,H()*1) operator* (const SCALE_MAT<TV>& a,const GRADIENT_VEC<TV,H>& u){return u*a.x;}

template<class TV,class H>
GRADIENT_VEC<TV,CM_FLAGS<H::n,H::mask|H::sp_mask,0> > operator* (const MATRIX<typename TV::SCALAR,TV::m>& a,const GRADIENT_VEC<TV,H>& u)
{
    GRADIENT_VEC<TV,CM_FLAGS<H::n,H::mask|H::sp_mask,0> > cm;
    Times_Helper(cm.c,u.c,a.Transposed());
    return cm;
}
}
}
#endif
