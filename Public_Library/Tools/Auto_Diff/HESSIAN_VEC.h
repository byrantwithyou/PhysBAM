//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HESSIAN_VEC
//##################################################################### 
#ifndef __HESSIAN_VEC__
#define __HESSIAN_VEC__

#include <Tools/Auto_Diff/GRADIENT_VEC.h>
#include <Tools/Auto_Diff/HESSIAN.h>
#include <Tools/Auto_Diff/MAT_TEN.h>
#include <Tools/Auto_Diff/PRIMITIVE_MATRICES.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{
template<class TV,class R> struct HESSIAN_VEC;

#define MK_HV(TV,FLAGS) HESSIAN_VEC<TV,decltype(FLAGS)>

template<class TV,class R>
struct HESSIAN_VEC
{
    typedef typename TV::SCALAR T;
    MK_MT(TV,R()) h;

    HESSIAN_VEC operator+ () const {return *this;}

    MK_HV(TV,-R()) operator- () const
    {MK_HV(TV,-R()) r;r.h.Neg(h);return r;}

    MK_HV(TV,R()*1) operator* (T a) const
    {MK_HV(TV,R()*1) r;r.h.Scale(h,a);return r;}

    MK_HV(TV,R()/1) operator/ (T a) const
    {return *this*(1/a);}

    template<class R2> MK_HV(TV,R()+R2())
    operator+ (const HESSIAN_VEC<TV,R2>& z) const
    {MK_HV(TV,R()+R2()) r;r.h.Add(h,z.h);return r;}

    template<class R2> MK_HV(TV,R()-R2())
    operator- (const HESSIAN_VEC<TV,R2>& z) const
    {MK_HV(TV,R()-R2()) r;r.h.Sub(h,z.h);return r;}

    template<class R2>
    void Fill_From(const HESSIAN_VEC<TV,R2>& z)
    {h.Fill_From(z.h);}
};

template<class TV,class R>
MK_HV(TV,R()*1) operator* (typename TV::SCALAR a,const HESSIAN_VEC<TV,R>& h)
{return h*a;}

template<class TV,class H>
MK_HV(TV,Tp0(V_FLAGS<1>(),H())) Tensor_Product_0(const HESSIAN<TV,H>& m,const TV& v)
{
    MK_HV(TV,Tp0(V_FLAGS<1>(),H())) hv;
    hv.h.Tensor_Product_0(m.h,v);
    return hv;
}

template<class TV,class H>
HESSIAN_VEC<TV,MT_FLAGS<H::n,0,0,0> > Tensor_Product_0(const HESSIAN<TV,H>& m,ZERO_VEC<TV> v)
{return HESSIAN_VEC<TV,MT_FLAGS<H::n,0,0,0> >();}

template<class TV>
void Symmetric_Transpose_Tensor_Product_12(MAT_TEN<TV,MT_EMPTY>& mt,const COL_MAT<TV,CM_EMPTY>& cm,const COL_VEC<TV,CV_EMPTY>& gv) {}

template<class TV,class G,class H,class R>
void Symmetric_Transpose_Tensor_Product_12(MAT_TEN<TV,R>& mt,const COL_MAT<TV,H>& cm,const COL_VEC<TV,G>& gv)
{
    mt.x=Symmetric_Tensor_Product_12(cm.x.Transposed(),gv.x);

    MK_CT(TV,Tp1(G::Base(),H::First())) c1;
    MK_CT(TV,Tp2(G::First(),H::Base())) c2;
    c1.Tensor_Product_1(cm.x.Transposed(),gv.Base());
    c2.Transpose_Tensor_Product_2(cm.Base(),gv.x);
    mt.c.Add(c1,c2);

    Symmetric_Transpose_Tensor_Product_12(mt.Base(),cm.Base(),gv.Base());
}

template<class TV,class G,class H>
MK_HV(TV,Stp12(G(),H())) Symmetric_Tensor_Product_12(const GRADIENT_VEC<TV,H>& gv,const GRADIENT<TV,G>& g)
{
    MK_HV(TV,Stp12(G(),H())) hv;
    Symmetric_Transpose_Tensor_Product_12(hv.h,gv.c,g.c);
    return hv;
}

template<class TV,class R> HESSIAN<TV,MM_FLAGS<R::n,0,0> > Contract_0(const HESSIAN_VEC<TV,R>& h,const ZERO_VEC<TV>& v) {return HESSIAN<TV,MM_FLAGS<R::n,0,0> >();}

template<class TV,class R> MK_H(TV,C0(R(),V_FLAGS<1>()))
Contract_0(const HESSIAN_VEC<TV,R>& h,const TV& v)
{
    MK_H(TV,C0(R(),V_FLAGS<1>())) o;
    Contract_0(o.h,h.h,v);
    return o;
}

template<class T_TENS,class H1,class H2> struct SYMMETRIC_DOUBLE_CONTRACT_HELPER;

template<> struct SYMMETRIC_DOUBLE_CONTRACT_HELPER<T_FLAGS<0,1,1>,CM_EMPTY,CM_EMPTY> {typedef MT_EMPTY TYPE;};

template<class H1,class H2>
struct SYMMETRIC_DOUBLE_CONTRACT_HELPER<T_FLAGS<0,1,1>,H1,H2>
{
    typedef decltype(Sdc12(T_FLAGS<0,1,1>(),H1::First(),H2::First())) X;
    typedef decltype(C2(T_FLAGS<0,1,1>(),H1::First())) CT1;
    typedef decltype(C2(T_FLAGS<0,1,1>(),H2::First())) CT2;
    typedef decltype(C1t(CT2(),H1::Base())+C1t(CT1(),H2::Base())) C;
    typedef typename SYMMETRIC_DOUBLE_CONTRACT_HELPER<T_FLAGS<0,1,1>,decltype(H1::Base()),decltype(H2::Base())>::TYPE B;
    typedef decltype(B().Append(X(),C())) TYPE;
};

template<class TV,class H,class H2> MK_HV(TV,(typename SYMMETRIC_DOUBLE_CONTRACT_HELPER<T_FLAGS<0,1,1>,H,H2>::TYPE()))
Symmetric_Double_Contract_12_With_Tensor(const PERM_TENSOR<TV>& t,const GRADIENT_VEC<TV,H>& a,const GRADIENT_VEC<TV,H2>& b)
{
    MK_HV(TV,(typename SYMMETRIC_DOUBLE_CONTRACT_HELPER<T_FLAGS<0,1,1>,H,H2>::TYPE())) hv;
    Symmetric_Double_Contract_12_With_Transposes(hv.h,t,a.c,b.c);
    return hv;
}

template<class TV,class R> HESSIAN_VEC<TV,MT_FLAGS<R::n,0,0,0> >
Contract_0(const HESSIAN_VEC<TV,R>& h,const ZERO_MAT<TV>& m)
{return HESSIAN_VEC<TV,MT_FLAGS<R::n,0,0,0> >();}

template<class TV,class R> HESSIAN_VEC<TV,R>
Contract_0(const HESSIAN_VEC<TV,R>& h,const ID_MAT<TV>& m)
{return h;}

template<class TV,class R> HESSIAN_VEC<TV,R>
Contract_0(const HESSIAN_VEC<TV,R>& h,const SCALE_MAT<TV>& m)
{return h*m.x;}

template<class TV,class R> MK_HV(TV,C0(R(),M_FLAGS<1,0>()))
Contract_0(const HESSIAN_VEC<TV,R>& h,const MATRIX<typename TV::SCALAR,TV::m>& m)
{
    MK_HV(TV,C0(R(),M_FLAGS<1,0>())) o;
    Contract_0(o.h,h.h,m);
    return o;
}

template<class R> MT_FLAGS<R::n,0,0,0> C0(R,M_FLAGS<0,0>);
template<class R> R C0(R,M_FLAGS<0,1>);
template<class R> R C0(R,M_FLAGS<1,1>);
}
}


#endif
