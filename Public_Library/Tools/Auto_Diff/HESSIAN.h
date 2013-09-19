//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HESSIAN
//##################################################################### 
#ifndef __HESSIAN__
#define __HESSIAN__

#include <Tools/Auto_Diff/GRADIENT.h>
#include <Tools/Auto_Diff/GRADIENT_VEC.h>
#include <Tools/Auto_Diff/MAT_MAT.h>
#include <Tools/Auto_Diff/PRIMITIVE_MATRICES.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{
template<class TV,class H> struct HESSIAN;

#define MK_H(TV,FLAGS) HESSIAN<TV,decltype(FLAGS)>

template<class TV,class H>
struct HESSIAN
{
    typedef typename TV::SCALAR T;
    MK_MM(TV,H()) h;

    HESSIAN operator+ () const {return *this;}

    MK_H(TV,-H()) operator- () const
    {MK_H(TV,-H()) r;r.h.Neg(h);return r;}

    MK_H(TV,H()*1) operator* (T a) const
    {MK_H(TV,H()*1) r;r.h.Scale(h,a);return r;}

    MK_H(TV,H()/1) operator/ (T a) const
    {return *this*(1/a);}

    template<class H2> MK_H(TV,H()+H2())
    operator+ (const HESSIAN<TV,H2>& z) const
    {MK_H(TV,H()+H2()) r;r.h.Add(h,z.h);return r;}

    template<class H2> MK_H(TV,H()-H2())
    operator- (const HESSIAN<TV,H2>& z) const
    {MK_H(TV,H()-H2()) r;r.h.Sub(h,z.h);return r;}

    template<class H2>
    void Fill_From(const HESSIAN<TV,H2>& z)
    {h.Fill_From(z.h);}

    MATRIX<T,TV::m> operator()(int i,int j) const {if(i>=j) return h(i,j);return h(j,i).Transposed();}

    SYMMETRIC_MATRIX<T,TV::m> operator()(int i) const {return h(i);}
};

template<class TV,class H>
MK_H(TV,H()*1) operator* (typename TV::SCALAR a,const HESSIAN<TV,H>& h)
{return h*a;}

template<class TV,class G,class G2> MK_H(TV,Sym_Outer(G(),G2()))
Symmetric_Outer_Product(const GRADIENT<TV,G>& u,const GRADIENT<TV,G2>& v)
{
    MK_H(TV,Sym_Outer(G(),G2())) H;
    Symmetric_Outer_Product_Helper(H.h,u.c,v.c);
    return H;
}

template<class TV,class G> MK_H(TV,Outer(G()))
Outer_Product(const GRADIENT<TV,G>& u)
{
    MK_H(TV,Outer(G())) H;
    Outer_Product_Helper(H.h,u.c);
    return H;
}

template<class TV,class H,class H2> MK_H(TV,Stt(H(),H2()))
Symmetric_Transpose_Times(const GRADIENT_VEC<TV,H>& u,const GRADIENT_VEC<TV,H2>& v)
{
    MK_H(TV,Stt(H(),H2())) h;
    Symmetric_Times_Transpose(h.h,u.c,v.c);
    return h;
}

template<class TV,class H> MK_H(TV,Tst(H())) Transpose_Times_Self(const GRADIENT_VEC<TV,H>& u)
{
    MK_H(TV,Tst(H())) h;
    Times_Self_Transpose(h.h,u.c);
    return h;
}

}
}


#endif
