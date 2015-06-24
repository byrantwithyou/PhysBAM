//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AUTO_HESS_EXT
//##################################################################### 
#ifndef __AUTO_HESS_EXT__
#define __AUTO_HESS_EXT__

#include <Tools/Auto_Diff/GRADIENT_MAT.h>
#include <Tools/Auto_Diff/HESSIAN_MAT.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/DIAGONAL_TENSOR.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{
using ::std::sin;
using ::std::cos;
using ::std::tan;
using ::std::exp;
using ::std::log;
using ::std::abs;
using ::std::sqrt;
using ::PhysBAM::sqr;
using ::PhysBAM::cube;

struct DIFF_UNUSED
{
    DIFF_UNUSED operator- () const {return DIFF_UNUSED();}
};

template<class T> struct DIFF_UNUSED_OK {static const int value=0;};
template<class T,class VEC> struct DIFF_UNUSED_OK<GRADIENT<T,VEC> > {static const int value=2;};
template<class T,class VEC> struct DIFF_UNUSED_OK<GRADIENT_VEC<T,VEC> > {static const int value=2;};
template<class T,class MAT> struct DIFF_UNUSED_OK<HESSIAN<T,MAT> > {static const int value=2;};
template<class T,class MAT> struct DIFF_UNUSED_OK<HESSIAN_VEC<T,MAT> > {static const int value=2;};
template<> struct DIFF_UNUSED_OK<DIFF_UNUSED> {static const int value=1;};

template<class T,class U>
inline typename enable_if<(DIFF_UNUSED_OK<T>::value | DIFF_UNUSED_OK<U>::value)==1,DIFF_UNUSED>::type
operator+ (const T&,const U&) {return DIFF_UNUSED();}

template<class T,class U>
inline typename enable_if<(DIFF_UNUSED_OK<T>::value | DIFF_UNUSED_OK<U>::value)==1,DIFF_UNUSED>::type
operator- (const T&,const U&) {return DIFF_UNUSED();}

template<class T,class U>
inline typename enable_if<(DIFF_UNUSED_OK<T>::value | DIFF_UNUSED_OK<U>::value)==1,DIFF_UNUSED>::type
operator* (const T&,const U&) {return DIFF_UNUSED();}

template<class T,class U>
inline typename enable_if<(DIFF_UNUSED_OK<T>::value | DIFF_UNUSED_OK<U>::value)==1,DIFF_UNUSED>::type
operator/ (const T&,const U&) {return DIFF_UNUSED();}

template<class T,class U>
inline typename enable_if<(DIFF_UNUSED_OK<T>::value | DIFF_UNUSED_OK<U>::value)==1,DIFF_UNUSED>::type
Contract_00(const T&,const U&) {return DIFF_UNUSED();}
template<class T,class U>
inline typename enable_if<(DIFF_UNUSED_OK<T>::value | DIFF_UNUSED_OK<U>::value)==1,DIFF_UNUSED>::type
Contract_10(const T&,const U&) {return DIFF_UNUSED();}
template<class T,class U>
inline typename enable_if<(DIFF_UNUSED_OK<T>::value | DIFF_UNUSED_OK<U>::value)==1,DIFF_UNUSED>::type
Contract_11(const T&,const U&) {return DIFF_UNUSED();}
template<class T,class U>
inline typename enable_if<(DIFF_UNUSED_OK<T>::value | DIFF_UNUSED_OK<U>::value)==1,DIFF_UNUSED>::type
Tensor_Product_01(const T&,const U&) {return DIFF_UNUSED();}
template<class T,class U>
inline typename enable_if<(DIFF_UNUSED_OK<T>::value | DIFF_UNUSED_OK<U>::value)==1,DIFF_UNUSED>::type
Double_Contract_00_11(const T&,const U&) {return DIFF_UNUSED();}
template<class T> inline DIFF_UNUSED Contract_0(const DIFF_UNUSED&,const T&) {return DIFF_UNUSED();}
template<class T> inline DIFF_UNUSED Tensor_Product_0(const DIFF_UNUSED&,const T&) {return DIFF_UNUSED();}
inline DIFF_UNUSED Choose(const DIFF_UNUSED&,const DIFF_UNUSED&) {return DIFF_UNUSED();}
inline void Fill_From(DIFF_UNUSED&,const DIFF_UNUSED&) {}

template<int Q,class T,class U> auto Symmetric_Outer_Product_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Symmetric_Outer_Product_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Symmetric_Outer_Product(t,u))>::type {return Symmetric_Outer_Product(t,u);}

template<int Q,class T,class U> auto Symmetric_Transpose_Times_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Symmetric_Transpose_Times_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Symmetric_Transpose_Times(t,u))>::type {return Symmetric_Transpose_Times(t,u);}

template<int Q,class T,class U> auto Outer_Product_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Outer_Product_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Outer_Product(t,u))>::type {return Outer_Product(t,u);}

template<int Q,class T,class U> auto Symmetric_Tensor_Product_12_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Symmetric_Tensor_Product_12_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Symmetric_Tensor_Product_12(t,u))>::type {return Symmetric_Tensor_Product_12(t,u);}

template<int Q,class T,class U> auto Symmetric_Tensor_Product_23_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Symmetric_Tensor_Product_23_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Symmetric_Tensor_Product_23(t,u))>::type {return Symmetric_Tensor_Product_23(t,u);}

template<int Q,class T,class U> auto Tensor_Product_2_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Tensor_Product_2_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Tensor_Product_2(t,u))>::type {return Tensor_Product_2(t,u);}

template<int Q,class T,class U> auto Tensor_Product_01_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Tensor_Product_01_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Tensor_Product_01(t,u))>::type {return Tensor_Product_01(t,u);}

template<int Q,class T,class U> auto Symmetric_Contract_10_12_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Symmetric_Contract_10_12_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Symmetric_Contract_10_12(t,u))>::type {return Symmetric_Contract_10_12(t,u);}

template<int Q,class T,class U> auto Symmetric_Contract_00_12_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Symmetric_Contract_00_12_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Symmetric_Contract_00_12(t,u))>::type {return Symmetric_Contract_00_12(t,u);}

template<int Q,class T,class U> auto Tensor_Product_1_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Tensor_Product_1_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Tensor_Product_1(t,u))>::type {return Tensor_Product_1(t,u);}

template<int Q,class T,class U> auto Symmetric_Contract_00_23_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Symmetric_Contract_00_23_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Symmetric_Contract_00_23(t,u))>::type {return Symmetric_Contract_00_23(t,u);}

template<int Q,class T,class U> auto Contract_01_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Contract_01_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Contract_01(t,u))>::type {return Contract_01(t,u);}

template<int Q,class T,class U> auto Symmetric_Contract_11_23_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Symmetric_Contract_11_23_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Symmetric_Contract_11_23(t,u))>::type {return Symmetric_Contract_11_23(t,u);}

template<int Q,class T,class U> auto Symmetric_Contract_01_23_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Symmetric_Contract_01_23_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Symmetric_Contract_01_23(t,u))>::type {return Symmetric_Contract_01_23(t,u);}

template<int Q,class T,class U> auto Symmetric_Double_Contract_00_11_01_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Symmetric_Double_Contract_00_11_01_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Symmetric_Double_Contract_00_11_01(t,u))>::type {return Symmetric_Double_Contract_00_11_01(t,u);}

template<int Q,class T> auto Contract_01_Q(const T&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T> auto Contract_01_Q(const T&t) -> typename enable_if<Q,decltype(Contract_01(t))>::type {return Contract_01_Q<Q&1>(t);}

template<int Q,class T> auto Transposed_01_Q(const T&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T> auto Transposed_01_Q(const T&t) -> typename enable_if<Q,decltype(Transposed_01(t))>::type {return Transposed_01_Q<Q&1>(t);}

template<int Q,class T> auto Twice_Symmetric_Part_23_Q(const T&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T> auto Twice_Symmetric_Part_23_Q(const T&t) -> typename enable_if<Q,decltype(Twice_Symmetric_Part_23(t))>::type {return Twice_Symmetric_Part_23(t);}

template<int Q,class T> auto Twice_Symmetric_Part_01_Q(const T&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T> auto Twice_Symmetric_Part_01_Q(const T&t) -> typename enable_if<Q,decltype(Twice_Symmetric_Part_01(t))>::type {return Twice_Symmetric_Part_01(t);}

template<int Q,class T,class U> auto Symmetric_Tensor_Product_02_23_Q(const T&,const U&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U> auto Symmetric_Tensor_Product_02_23_Q(const T&t,const U&u) -> typename enable_if<Q,decltype(Symmetric_Tensor_Product_02_23(t,u))>::type {return Symmetric_Tensor_Product_02_23(t,u);}

template<int Q,class T> auto Transpose_Times_Self_Q(const T&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T> auto Transpose_Times_Self_Q(const T&t) -> typename enable_if<Q,decltype(Transpose_Times_Self(t))>::type {return Transpose_Times_Self(t);}

template<int Q,class T> auto Outer_Product_Q(const T&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T> auto Outer_Product_Q(const T&t) -> typename enable_if<Q,decltype(Outer_Product(t))>::type {return Outer_Product(t);}

template<int Q,class T,class U,class V> auto Symmetric_Double_Contract_12_With_Tensor_Q(const T&,const U&,const V&) -> typename enable_if<!(Q),DIFF_UNUSED>::type {return DIFF_UNUSED();}
template<int Q,class T,class U,class V> auto Symmetric_Double_Contract_12_With_Tensor_Q(const T&t,const U&u,const V&v) -> typename enable_if<Q,decltype(Symmetric_Double_Contract_12_With_Tensor(t,u,v))>::type {return Symmetric_Double_Contract_12_With_Tensor(t,u,v);}



template<class T,class VEC,class MAT,int Q> struct AUTO_HESS_EXT;

template<class T,class VEC,class MAT> AUTO_HESS_EXT<T,VEC,MAT,3>
Make_Hess(T x,const GRADIENT<T,VEC>& dx,const HESSIAN<T,MAT>& ddx);

template<class T,class VEC> AUTO_HESS_EXT<T,VEC,DIFF_UNUSED,1>
Make_Hess(T x,const GRADIENT<T,VEC>& dx,const DIFF_UNUSED& ddx);

template<class T> AUTO_HESS_EXT<T,DIFF_UNUSED,DIFF_UNUSED,0>
Make_Hess(T x,const DIFF_UNUSED& dx,const DIFF_UNUSED& ddx);

template<class T,class VEC,class MAT,int Q>
struct AUTO_HESS_EXT
{
    static_assert(is_scalar<T>::value,"This definition of AUTO_HESS_EXT is only for scalars");
    typedef typename conditional<Q&1,GRADIENT<T,VEC>,DIFF_UNUSED>::type GRAD;
    typedef typename conditional<Q&2,HESSIAN<T,MAT>,DIFF_UNUSED>::type HESS;

    T x;
    GRAD dx;
    HESS ddx;

    AUTO_HESS_EXT(T z=T()): x(z)
    {}

    AUTO_HESS_EXT(T z,const GRAD& dz,const HESS& ddz)
        :x(z),dx(dz),ddx(ddz)
    {}

    decltype(Make_Hess(-x,-dx,-ddx)) operator-() const
    {return Make_Hess(-x,-dx,-ddx);}

    AUTO_HESS_EXT operator+() const
    {return *this;}

    template<class VEC1,class MAT1> auto
    operator+(const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess(this->x+a.x,this->dx+a.dx,this->ddx+a.ddx))
    {return Make_Hess(x+a.x,dx+a.dx,ddx+a.ddx);}

    template<class VEC1,class MAT1> auto
    operator-(const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess(x-a.x,dx-a.dx,ddx-a.ddx))
    {return Make_Hess(x-a.x,dx-a.dx,ddx-a.ddx);}

    template<class VEC1,class MAT1> auto
    operator*(const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess(this->x*a.x,a.x*this->dx+this->x*a.dx,a.x*this->ddx+this->x*a.ddx+Symmetric_Outer_Product_Q<Q&2>(this->dx,a.dx)))
    {return Make_Hess(x*a.x,a.x*dx+x*a.dx,a.x*ddx+x*a.ddx+Symmetric_Outer_Product_Q<Q&2>(dx,a.dx));}

    template<class VEC1,class MAT1> auto
    operator/(const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess(T(),(this->dx-T()*a.dx)/a.x,this->ddx/a.x-T()*a.ddx/a.x-Symmetric_Outer_Product_Q<Q&2>((this->dx-T()*a.dx),a.dx/T())))
    {
        T z=x/a.x,ax2=sqr(a.x);
        auto q=dx-z*a.dx;
        return Make_Hess(z,q/a.x,ddx/a.x-z*a.ddx/a.x-Symmetric_Outer_Product_Q<Q&2>(q,a.dx/ax2));
    }

    AUTO_HESS_EXT operator+(T a) const
    {return Make_Hess(x+a,dx,ddx);}

    AUTO_HESS_EXT operator-(T a) const
    {return Make_Hess(x-a,dx,ddx);}

    auto operator*(T a) const -> decltype(Make_Hess(this->x/a,this->dx/a,this->ddx/a))
    {return Make_Hess(x*a,a*dx,a*ddx);}

    auto operator/(T a) const -> decltype(Make_Hess(this->x/a,this->dx/a,this->ddx/a))
    {return Make_Hess(x/a,dx/a,ddx/a);}

    template<class VEC1,class MAT1>
    void Fill_From(const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& z)
    {x=z.x;::PhysBAM::HETERO_DIFF::Fill_From(dx,z.dx);::PhysBAM::HETERO_DIFF::Fill_From(ddx,z.ddx);}
};

template<class T,class VEC,class MAT> AUTO_HESS_EXT<T,VEC,MAT,3>
Make_Hess(T x,const GRADIENT<T,VEC>& dx,const HESSIAN<T,MAT>& ddx)
{return AUTO_HESS_EXT<T,VEC,MAT,3>(x,dx,ddx);}

template<class T,class VEC> AUTO_HESS_EXT<T,VEC,DIFF_UNUSED,1>
Make_Hess(T x,const GRADIENT<T,VEC>& dx,const DIFF_UNUSED& ddx)
{return AUTO_HESS_EXT<T,VEC,DIFF_UNUSED,1>(x,dx,ddx);}

template<class T> AUTO_HESS_EXT<T,DIFF_UNUSED,DIFF_UNUSED,0>
Make_Hess(T x,const DIFF_UNUSED& dx,const DIFF_UNUSED& ddx)
{return AUTO_HESS_EXT<T,DIFF_UNUSED,DIFF_UNUSED,0>(x,dx,ddx);}

template<class LAYOUT,class T>
typename enable_if<is_scalar<T>::value,AUTO_HESS_EXT<T,typename EMPTY_VEC<1,LAYOUT,-1>::TYPE,typename EMPTY_MAT<LAYOUT,-1>::TYPE,3> >::type
Hess_From_Const(T a)
{return AUTO_HESS_EXT<T,typename EMPTY_VEC<1,LAYOUT,-1>::TYPE,typename EMPTY_MAT<LAYOUT,-1>::TYPE,3>(a);}

template<class LAYOUT,class T>
typename enable_if<is_scalar<T>::value,AUTO_HESS_EXT<T,typename EMPTY_VEC<1,LAYOUT,-1>::TYPE,DIFF_UNUSED,1> >::type
Diff_From_Const(T a)
{return AUTO_HESS_EXT<T,typename EMPTY_VEC<1,LAYOUT,-1>::TYPE,DIFF_UNUSED,1>(a);}

template<class LAYOUT,class T>
typename enable_if<is_scalar<T>::value,AUTO_HESS_EXT<T,DIFF_UNUSED,DIFF_UNUSED,0> >::type
From_Const(T a)
{return AUTO_HESS_EXT<T,DIFF_UNUSED,DIFF_UNUSED,0>(a);}

template<class T,class VEC,class MAT,int Q>
inline AUTO_HESS_EXT<T,VEC,MAT,Q> operator+(T c,const AUTO_HESS_EXT<T,VEC,MAT,Q>& a)
{return Make_Hess(c+a.x,a.dx,a.ddx);}

template<class T,class VEC,class MAT,int Q>
inline auto operator-(T c,const AUTO_HESS_EXT<T,VEC,MAT,Q>& a) -> decltype(Make_Hess(c-a.x,-a.dx,-a.ddx))
{return Make_Hess(c-a.x,-a.dx,-a.ddx);}

template<class T,class VEC,class MAT,int Q> inline auto
operator*(T c,const AUTO_HESS_EXT<T,VEC,MAT,Q>& a) -> decltype(Make_Hess(c*a.x,c*a.dx,c*a.ddx))
{return Make_Hess(c*a.x,c*a.dx,c*a.ddx);}

template<class T,class VEC,class MAT,int Q> inline auto
operator/(T c,const AUTO_HESS_EXT<T,VEC,MAT,Q>& a)
    -> decltype(Make_Hess(T(),-T()*a.dx,2*T()/a.x*Outer_Product_Q<Q&2>(a.dx)-T()*a.ddx))
{T z=c/a.x,w=z/a.x;return Make_Hess(z,-w*a.dx,2*w/a.x*Outer_Product_Q<Q&2>(a.dx)-w*a.ddx);}

template<class T,class VEC,class MAT,int Q> inline auto
sqrt(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a)
    -> decltype(Make_Hess(T(),(a.dx/(2*T())),a.ddx/(2*T())-Outer_Product_Q<Q&2>(a.dx/(2*T()))/T()))
{T s=sqrt(a.x);auto t=a.dx/(2*s);return Make_Hess(s,t,a.ddx/(2*s)-Outer_Product_Q<Q&2>(t)/s);}

template<class T,class VEC,class MAT,int Q> inline auto
sqr(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a)
    -> decltype(Make_Hess(sqr(a.x),2*a.x*a.dx,2*a.x*a.ddx+Outer_Product_Q<Q&2>(a.dx)*2))
{return Make_Hess(sqr(a.x),2*a.x*a.dx,2*a.x*a.ddx+Outer_Product_Q<Q&2>(a.dx)*2);}

template<class T,class VEC,class MAT,int Q> inline auto
cube(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a)
    -> decltype(Make_Hess(a.x*T(),3*T()*a.dx,3*T()*a.ddx+Outer_Product_Q<Q&2>(a.dx)*(6*a.x)))
{T sq=sqr(a.x);return Make_Hess(a.x*sq,3*sq*a.dx,3*sq*a.ddx+Outer_Product_Q<Q&2>(a.dx)*(6*a.x));}

template<class T,class VEC,class MAT,int Q,class VEC1,class MAT1> inline
AUTO_HESS_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC1())),decltype(MAT_CHOOSE::Type(MAT(),MAT1())),Q>
max(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& b)
{
    AUTO_HESS_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC1())),decltype(MAT_CHOOSE::Type(MAT(),MAT1())),Q> r;
    if(a.x>b.x) r.Fill_From(a);
    else r.Fill_From(b);
    return r;
}

template<class T,class VEC,class MAT,int Q> inline auto
max(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,T b) -> AUTO_HESS_EXT<T,VEC,MAT,Q>
{
    if(a.x>b) return a;
    AUTO_HESS_EXT<T,VEC,MAT,Q> r;
    r.x=b;
    return r;
}

template<class T,class VEC,class MAT,int Q> inline auto
max(T b,const AUTO_HESS_EXT<T,VEC,MAT,Q>& a) -> AUTO_HESS_EXT<T,VEC,MAT,Q>
{
    if(a.x>b) return a;
    AUTO_HESS_EXT<T,VEC,MAT,Q> r;
    r.x=b;
    return r;
}

template<class T,class VEC,class MAT,int Q>
inline AUTO_HESS_EXT<T,VEC,MAT,Q>
max(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,const AUTO_HESS_EXT<T,VEC,MAT,Q>& b)
{return a.x>b.x?a:b;}

template<class T,class VEC,class MAT,int Q,class VEC1,class MAT1> inline
AUTO_HESS_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC1())),decltype(MAT_CHOOSE::Type(MAT(),MAT1())),Q>
min(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& b)
{
    AUTO_HESS_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC1())),decltype(MAT_CHOOSE::Type(MAT(),MAT1())),Q> r;
    if(a.x>b.x) r.Fill_From(b);
    else r.Fill_From(a);
    return r;
}

template<class T,class VEC,class MAT,int Q> inline auto
min(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,T b) -> AUTO_HESS_EXT<T,VEC,MAT,Q>
{
    if(a.x<b) return a;
    AUTO_HESS_EXT<T,VEC,MAT,Q> r;
    r.x=b;
    return r;
}

template<class T,class VEC,class MAT,int Q> inline auto
min(T b,const AUTO_HESS_EXT<T,VEC,MAT,Q>& a) -> AUTO_HESS_EXT<T,VEC,MAT,Q>
{
    if(a.x<b) return a;
    AUTO_HESS_EXT<T,VEC,MAT,Q> r;
    r.x=b;
    return r;
}

template<class T,class VEC,class MAT,int Q>
inline AUTO_HESS_EXT<T,VEC,MAT,Q>
min(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,const AUTO_HESS_EXT<T,VEC,MAT,Q>& b)
{return a.x>b.x?b:a;}

template<class T,class VEC,class MAT,int Q> inline auto
log(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a)
    -> decltype(Make_Hess(::std::log(a.x),a.dx/a.x,a.ddx/a.x-Outer_Product_Q<Q&2>(a.dx/a.x)))
{auto z=a.dx/a.x;return Make_Hess(::std::log(a.x),z,a.ddx/a.x-Outer_Product_Q<Q&2>(z));}

template<class T,class VEC,class MAT,int Q> inline auto
exp(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a)
    -> decltype(Make_Hess(T(),T()*a.dx,T()*(a.ddx+Outer_Product_Q<Q&2>(a.dx))))
{T s=exp(a.x);return Make_Hess(s,s*a.dx,s*(a.ddx+Outer_Product_Q<Q&2>(a.dx)));}

template<class T,class VEC,class MAT,int Q> inline auto
sin(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a)
    -> decltype(Make_Hess(T(),T()*a.dx,T()*a.ddx-T()*Outer_Product_Q<Q&2>(a.dx)))
{T s=sin(a.x),c=cos(a.x);return Make_Hess(s,c*a.dx,c*a.ddx-s*Outer_Product_Q<Q&2>(a.dx));}

template<class T,class VEC,class MAT,int Q> inline auto
cos(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a)
    -> decltype(Make_Hess(T(),-T()*a.dx,-T()*a.ddx-T()*Outer_Product_Q<Q&2>(a.dx)))
{T s=sin(a.x),c=cos(a.x);return Make_Hess(c,-s*a.dx,-s*a.ddx-c*Outer_Product_Q<Q&2>(a.dx));}

template<class T,class VEC,class MAT,int Q> inline auto
tan(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a)
    -> decltype(Make_Hess(T(),T()*a.dx,T()*a.ddx+2*T()*T()*Outer_Product_Q<Q&2>(a.dx)))
{T t=tan(a.x),s=1+t*t;return Make_Hess(t,s*a.dx,s*a.ddx+2*s*t*Outer_Product_Q<Q&2>(a.dx));}

template<class T,class VEC,class MAT,int Q> inline
AUTO_HESS_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC_NEG::Type(VEC()))),decltype(MAT_CHOOSE::Type(MAT(),MAT_NEG::Type(MAT()))),Q>
abs(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a)
{
    AUTO_HESS_EXT<T,decltype(VEC_CHOOSE::Type(VEC(),VEC_NEG::Type(VEC()))),decltype(MAT_CHOOSE::Type(MAT(),MAT_NEG::Type(MAT()))),Q> r;
    if(a.x>=0) r.Fill_From(a);
    else r.Fill_From(-a);
    return r;
}

template<class T,class VEC,class MAT,int Q,class VEC1,class MAT1> inline auto
hypot(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& b)
    -> decltype(Make_Hess(T(),T()*a.dx+T()*b.dx,T()*a.ddx+T()*b.ddx+Outer_Product_Q<Q&2>(T()*b.dx-T()*a.dx)/T()))
{
    T c=sqrt(sqr(a.x)+sqr(b.x)),d=a.x/c,e=b.x/c;
    return Make_Hess(c,d*a.dx+e*b.dx,d*a.ddx+e*b.ddx+Outer_Product_Q<Q&2>(d*b.dx-e*a.dx)/c);
}

template<class T,class VEC,class MAT,int Q> inline auto
hypot(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,T b)
    -> decltype(Make_Hess(T(),T()*a.dx,T()*a.ddx+Outer_Product_Q<Q&2>(T()*a.dx)/T()))
{
    T c=sqrt(sqr(a.x)+sqr(b)),d=a.x/c,e=b/c;
    return Make_Hess(c,d*a.dx,d*a.ddx+Outer_Product_Q<Q&2>(e*a.dx)/c);
}

template<class T,class VEC,class MAT,int Q> inline auto
hypot(T b,const AUTO_HESS_EXT<T,VEC,MAT,Q>& a)
    -> decltype(hypot(a,b))
{return hypot(a,b);}

template<class T,class VEC,class MAT,int Q,class VEC1,class MAT1,class VEC2,class MAT2> inline auto
hypot(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& b,const AUTO_HESS_EXT<T,VEC2,MAT2,Q>& c)
    -> decltype(Make_Hess(T(),T()*a.dx+T()*b.dx+T()*c.dx,T()*a.ddx+T()*b.ddx+T()*c.ddx+(Outer_Product_Q<Q&2>(T()*b.dx-T()*a.dx)+Outer_Product_Q<Q&2>(T()*c.dx-T()*b.dx)+Outer_Product_Q<Q&2>(T()*a.dx-T()*c.dx))/T()))
{
    T s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c.x)),aa=a.x/s,bb=b.x/s,cc=c.x/s;
    auto ab=Outer_Product_Q<Q&2>(aa*b.dx-bb*a.dx);
    auto bc=Outer_Product_Q<Q&2>(bb*c.dx-cc*b.dx);
    auto ca=Outer_Product_Q<Q&2>(cc*a.dx-aa*c.dx);
    return Make_Hess(s,aa*a.dx+bb*b.dx+cc*c.dx,aa*a.ddx+bb*b.ddx+cc*c.ddx+(ab+bc+ca)/s);
}

template<class T,class VEC,class MAT,int Q,class VEC1,class MAT1> inline auto
hypot(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& b,T c)
    -> decltype(Make_Hess(T(),T()*a.dx+T()*b.dx,T()*a.ddx+T()*b.ddx+Outer_Product_Q<Q&2>(T()*b.dx-T()*a.dx)/T()+(Outer_Product_Q<Q&2>(b.dx)+Outer_Product_Q<Q&2>(a.dx))*(sqr(T())/T())))
{
    T s=sqrt(sqr(a.x)+sqr(b.x)+sqr(c)),aa=a.x/s,bb=b.x/s,cc=c/s;
    auto ab=Outer_Product_Q<Q&2>(aa*b.dx-bb*a.dx);
    auto bc=Outer_Product_Q<Q&2>(b.dx);
    auto ca=Outer_Product_Q<Q&2>(a.dx);
    return Make_Hess(s,aa*a.dx+bb*b.dx,aa*a.ddx+bb*b.ddx+ab/s+(bc+ca)*(sqr(cc)/s));
}

template<class T,class VEC,class MAT,int Q,class VEC1,class MAT1> inline auto
hypot(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,T c,const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& b)
    -> decltype(hypot(a,b,c))
{return hypot(a,b,c);}

template<class T,class VEC,class MAT,int Q,class VEC1,class MAT1> inline auto
hypot(T c,const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& b)
    -> decltype(hypot(a,b,c))
{return hypot(a,b,c);}

template<class T,class VEC,class MAT,int Q> inline auto
hypot(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,T b,T c)
    -> decltype(Make_Hess(T(),T()*a.dx,T()*a.ddx+Outer_Product_Q<Q&2>(a.dx)*T()))
{
    T t=sqr(b)+sqr(c),s=sqrt(sqr(a.x)+t),d=a.x/s,e=t/(s*s*s);
    return Make_Hess(s,d*a.dx,d*a.ddx+Outer_Product_Q<Q&2>(a.dx)*e);
}

template<class T,class VEC,class MAT,int Q> inline auto
hypot(T b,const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,T c)
    -> decltype(hypot(a,b,c))
{return hypot(a,b,c);}

template<class T,class VEC,class MAT,int Q> inline auto
hypot(T b,T c,const AUTO_HESS_EXT<T,VEC,MAT,Q>& a)
    -> decltype(hypot(a,b,c))
{return hypot(a,b,c);}

template<class T,class VEC,class MAT,int Q,class VEC1,class MAT1> inline auto
atan2(const AUTO_HESS_EXT<T,VEC,MAT,Q>& y,const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& x)
    -> decltype(Make_Hess(::std::atan2(y.x,x.x),T()*y.dx-T()*x.dx,T()*y.ddx-T()*x.ddx-Symmetric_Outer_Product_Q<Q&2>(T()*y.dx-T()*x.dx,T()*x.dx+T()*y.dx)))
{
    T c=sqr(x.x)+sqr(y.x),d=x.x/c,e=y.x/c;
    auto f=d*y.dx-e*x.dx;
    return Make_Hess(::std::atan2(y.x,x.x),f,d*y.ddx-e*x.ddx-Symmetric_Outer_Product_Q<Q&2>(f,d*x.dx+e*y.dx));
}

template<class T,class VEC,class MAT,int Q> inline auto
atan2(const AUTO_HESS_EXT<T,VEC,MAT,Q>& y,T x)
    -> decltype(Make_Hess(::std::atan2(y.x,x),T()*y.dx,T()*y.ddx-Symmetric_Outer_Product_Q<Q&2>(T()*y.dx,T()*y.dx)))
{
    T c=sqr(x)+sqr(y.x),d=x/c,e=y.x/c;
    auto f=d*y.dx;
    return Make_Hess(::std::atan2(y.x,x),f,d*y.ddx-Symmetric_Outer_Product_Q<Q&2>(f,e*y.dx));
}

template<class T,class VEC,class MAT,int Q> inline auto
atan2(T y,const AUTO_HESS_EXT<T,VEC,MAT,Q>& x)
    -> decltype(Make_Hess(::std::atan2(y,x.x),-T()*x.dx,-T()*x.ddx-Symmetric_Outer_Product_Q<Q&2>(-T()*x.dx,T()*x.dx)))
{
    T c=sqr(x.x)+sqr(y),d=x.x/c,e=y/c;
    auto f=-e*x.dx;
    return Make_Hess(::std::atan2(y,x.x),f,-e*x.ddx-Symmetric_Outer_Product_Q<Q&2>(f,d*x.dx));
}
template<class TV,class VEC,class MAT,int Q> struct AUTO_HESS_EXT_VEC;

template<class TV,class VEC,class MAT> AUTO_HESS_EXT_VEC<TV,VEC,MAT,3>
Make_Hess_Vec(const TV& x,const GRADIENT_VEC<TV,VEC>& dx,const HESSIAN_VEC<TV,MAT>& ddx);

template<class TV,class VEC> AUTO_HESS_EXT_VEC<TV,VEC,DIFF_UNUSED,1>
Make_Hess_Vec(const TV& x,const GRADIENT_VEC<TV,VEC>& dx,const DIFF_UNUSED& ddx);

template<class TV> AUTO_HESS_EXT_VEC<TV,DIFF_UNUSED,DIFF_UNUSED,0>
Make_Hess_Vec(const TV& x,const DIFF_UNUSED& dx,const DIFF_UNUSED& ddx);

template<class T,class TV,class VEC1,class MAT1,int Q,class VEC2,class MAT2> auto
Cross_Helper(const AUTO_HESS_EXT_VEC<VECTOR<T,3>,VEC1,MAT1,Q>& a,const AUTO_HESS_EXT_VEC<TV,VEC2,MAT2,Q>& b)
    -> decltype(Make_Hess_Vec(TV().Cross(TV()),MATRIX<T,TV::m>()*typename conditional<Q&1,GRADIENT_VEC<TV,VEC2>,DIFF_UNUSED>::type()-MATRIX<T,TV::m>()*typename conditional<Q&1,GRADIENT_VEC<TV,VEC1>,DIFF_UNUSED>::type(),Contract_00(typename conditional<Q&2,HESSIAN_VEC<TV,MAT1>,DIFF_UNUSED>::type(),MATRIX<T,TV::m>())-Contract_00(typename conditional<Q&2,HESSIAN_VEC<TV,MAT2>,DIFF_UNUSED>::type(),MATRIX<T,TV::m>())+Symmetric_Double_Contract_12_With_Tensor_Q<Q&2>(PERMUTATION_TENSOR<T>(1),typename conditional<Q&1,GRADIENT_VEC<TV,VEC1>,DIFF_UNUSED>::type(),typename conditional<Q&1,GRADIENT_VEC<TV,VEC2>,DIFF_UNUSED>::type())));

template<class T,class VEC1,class MAT1,int Q> auto
Cross_Helper(const AUTO_HESS_EXT_VEC<VECTOR<T,3>,VEC1,MAT1,Q>& a,const VECTOR<T,3>& b)
    -> decltype(Make_Hess_Vec(VECTOR<T,3>(),-MATRIX<T,3>()*typename conditional<Q&1,GRADIENT_VEC<VECTOR<T,3>,VEC1>,DIFF_UNUSED>::type(),Contract_00(typename conditional<Q&2,HESSIAN_VEC<VECTOR<T,3>,MAT1>,DIFF_UNUSED>::type(),MATRIX<T,3>())));

template<class TV,class VEC1,class MAT1,int Q,class TYPE>
typename enable_if<TV::m!=3,AUTO_HESS_EXT_VEC<TV,VEC1,MAT1,Q> >::type
Cross_Helper(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1,Q>& a,TYPE);

template<class TV,class VEC,class MAT,int Q>
struct AUTO_HESS_EXT_VEC
{
    typedef typename TV::SCALAR T;
    typedef typename conditional<Q&1,GRADIENT_VEC<TV,VEC>,DIFF_UNUSED>::type GRAD;
    typedef typename conditional<Q&2,HESSIAN_VEC<TV,MAT>,DIFF_UNUSED>::type HESS;

    TV x;
    GRAD dx;
    HESS ddx;

    AUTO_HESS_EXT_VEC(TV x=TV()):
        x(x)
    {}

    AUTO_HESS_EXT_VEC(TV x,const GRAD& dx,const HESS& ddx):
        x(x),dx(dx),ddx(ddx)
    {}

    auto operator-() const -> decltype(Make_Hess_Vec(-this->x,-this->dx,-this->ddx))
    {return Make_Hess_Vec(-x,-dx,-ddx);}

    AUTO_HESS_EXT_VEC operator+() const
    {return *this;}

    template<class VEC1,class MAT1> auto
    operator+(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess_Vec(x+a.x,dx+a.dx,ddx+a.ddx))
    {return Make_Hess_Vec(x+a.x,dx+a.dx,ddx+a.ddx);}

    template<class VEC1,class MAT1> auto
    operator-(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess_Vec(x-a.x,dx-a.dx,ddx-a.ddx))
    {return Make_Hess_Vec(x-a.x,dx-a.dx,ddx-a.ddx);}

    template<class VEC1,class MAT1> auto
    operator*(const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess_Vec(this->x*a.x,a.x*this->dx+Outer_Product_Q<Q&1>(this->x,a.dx),this->ddx*a.x+Tensor_Product_0(a.ddx,this->x)+Symmetric_Tensor_Product_12_Q<Q&2>(this->dx,a.dx)))
    {return Make_Hess_Vec(x*a.x,a.x*dx+Outer_Product_Q<Q&1>(x,a.dx),ddx*a.x+Tensor_Product_0(a.ddx,x)+Symmetric_Tensor_Product_12_Q<Q&2>(dx,a.dx));}

    template<class VEC1,class MAT1> auto
    operator/(const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess_Vec(this->x/a.x,this->dx/a.x+Outer_Product_Q<Q&1>(x,-a.dx/(a.x*a.x)),this->ddx/a.x+Tensor_Product_0(2/a.x*Outer_Product_Q<Q&2>(a.dx)-a.ddx,this->x/(a.x*a.x))+Symmetric_Tensor_Product_12_Q<Q&2>(this->dx,-a.dx/(a.x*a.x))))
    {auto p=-a.dx/(a.x*a.x);return Make_Hess_Vec(x/a.x,dx/a.x+Outer_Product_Q<Q&1>(x,p),ddx/a.x+Tensor_Product_0(2/a.x*Outer_Product_Q<Q&2>(a.dx)-a.ddx,x/(a.x*a.x))+Symmetric_Tensor_Product_12_Q<Q&2>(dx,p));}

    AUTO_HESS_EXT_VEC operator+(TV a) const
    {return Make_Hess_Vec(x+a,dx,ddx);}

    AUTO_HESS_EXT_VEC operator-(TV a) const
    {return Make_Hess_Vec(x-a,dx,ddx);}

    auto operator*(T a) const -> decltype(Make_Hess_Vec(this->x*a,this->dx*a,this->ddx*a))
    {return Make_Hess_Vec(x*a,dx*a,ddx*a);}

    auto operator/(T a) const -> decltype(Make_Hess_Vec(this->x/a,this->dx/a,this->ddx/a))
    {return Make_Hess_Vec(x/a,dx/a,ddx/a);}

    auto Dot(TV v) const -> decltype(Make_Hess(this->x.Dot(v),this->dx.Transpose_Times(v),Contract_0(this->ddx,v)))
    {return Make_Hess(x.Dot(v),dx.Transpose_Times(v),Contract_0(ddx,v));}

    template<class VEC1,class MAT1> auto
    Dot(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess(this->x.Dot(a.x),this->dx.Transpose_Times(a.x)+a.dx.Transpose_Times(this->x),Contract_0(this->ddx,a.x)+Contract_0(a.ddx,this->x)+Symmetric_Transpose_Times_Q<Q&2>(this->dx,a.dx)))
    {return Make_Hess(x.Dot(a.x),dx.Transpose_Times(a.x)+a.dx.Transpose_Times(x),Contract_0(ddx,a.x)+Contract_0(a.ddx,x)+Symmetric_Transpose_Times_Q<Q&2>(dx,a.dx));}

    auto Cross(const TV& a) const -> decltype(Cross_Helper(*this,a))
    {return Cross_Helper(*this,a);}

    template<class VEC1,class MAT1>
    auto Cross(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1,Q>& a) const -> decltype(Cross_Helper(*this,a))
    {return Cross_Helper(*this,a);}

    decltype(Make_Hess(x.Magnitude_Squared(),dx.Transpose_Times(x)*2,(Contract_0(ddx,x)+Transpose_Times_Self_Q<Q&2>(dx))*2)) Magnitude_Squared() const
    {return Make_Hess(x.Magnitude_Squared(),dx.Transpose_Times(x)*2,(Contract_0(ddx,x)+Transpose_Times_Self_Q<Q&2>(dx))*2);}

    auto Magnitude() const
        -> decltype(Make_Hess(T(),dx.Transpose_Times(x)/T(),(Contract_0(ddx,x)+Transpose_Times_Self_Q<Q&2>(dx))/T()-Outer_Product_Q<Q&2>(dx.Transpose_Times(x)/T())/T()))
    {T s=x.Magnitude();auto t=dx.Transpose_Times(x)/s;return Make_Hess(s,t,(Contract_0(ddx,x)+Transpose_Times_Self_Q<Q&2>(dx))/s-Outer_Product_Q<Q&2>(t)/s);}
};


template<class T,class TV,class VEC1,class MAT1,int Q,class VEC2,class MAT2> auto
Cross_Helper(const AUTO_HESS_EXT_VEC<VECTOR<T,3>,VEC1,MAT1,Q>& a,const AUTO_HESS_EXT_VEC<TV,VEC2,MAT2,Q>& b)
    -> decltype(Make_Hess_Vec(TV().Cross(TV()),MATRIX<T,TV::m>()*typename conditional<Q&1,GRADIENT_VEC<TV,VEC2>,DIFF_UNUSED>::type()-MATRIX<T,TV::m>()*typename conditional<Q&1,GRADIENT_VEC<TV,VEC1>,DIFF_UNUSED>::type(),Contract_00(typename conditional<Q&2,HESSIAN_VEC<TV,MAT1>,DIFF_UNUSED>::type(),MATRIX<T,TV::m>())-Contract_00(typename conditional<Q&2,HESSIAN_VEC<TV,MAT2>,DIFF_UNUSED>::type(),MATRIX<T,TV::m>())+Symmetric_Double_Contract_12_With_Tensor_Q<Q&2>(PERMUTATION_TENSOR<T>(1),typename conditional<Q&1,GRADIENT_VEC<TV,VEC1>,DIFF_UNUSED>::type(),typename conditional<Q&1,GRADIENT_VEC<TV,VEC2>,DIFF_UNUSED>::type())))
{
    MATRIX<T,TV::m> cp_t=MATRIX<T,TV::m>::Cross_Product_Matrix(a.x),cp_a=MATRIX<T,TV::m>::Cross_Product_Matrix(b.x);
    return Make_Hess_Vec(a.x.Cross(b.x),cp_t*b.dx-cp_a*a.dx,Contract_00(a.ddx,cp_a)-Contract_00(b.ddx,cp_t)+Symmetric_Double_Contract_12_With_Tensor_Q<Q&2>(PERMUTATION_TENSOR<T>(1),a.dx,b.dx));
}

template<class T,class VEC1,class MAT1,int Q> auto
Cross_Helper(const AUTO_HESS_EXT_VEC<VECTOR<T,3>,VEC1,MAT1,Q>& a,const VECTOR<T,3>& b)
    -> decltype(Make_Hess_Vec(VECTOR<T,3>(),-MATRIX<T,3>()*typename conditional<Q&1,GRADIENT_VEC<VECTOR<T,3>,VEC1>,DIFF_UNUSED>::type(),Contract_00(typename conditional<Q&2,HESSIAN_VEC<VECTOR<T,3>,MAT1>,DIFF_UNUSED>::type(),MATRIX<T,3>())))
{
    MATRIX<T,3> cp_a=MATRIX<T,3>::Cross_Product_Matrix(b);
    return Make_Hess_Vec(a.x.Cross(b),-cp_a*a.dx,Contract_00(a.ddx,cp_a));
}

template<class TV,class VEC1,class MAT1,int Q,class TYPE>
typename enable_if<TV::m!=3,AUTO_HESS_EXT_VEC<TV,VEC1,MAT1,Q> >::type
Cross_Helper(const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1,Q>& a,TYPE)
{
    PHYSBAM_FATAL_ERROR("Cross product not defined except in 3D.");
    return a;
}

template<class TV,class VEC,class MAT> AUTO_HESS_EXT_VEC<TV,VEC,MAT,3>
Make_Hess_Vec(const TV& x,const GRADIENT_VEC<TV,VEC>& dx,const HESSIAN_VEC<TV,MAT>& ddx)
{return AUTO_HESS_EXT_VEC<TV,VEC,MAT,3>(x,dx,ddx);}

template<class TV,class VEC> AUTO_HESS_EXT_VEC<TV,VEC,DIFF_UNUSED,1>
Make_Hess_Vec(const TV& x,const GRADIENT_VEC<TV,VEC>& dx,const DIFF_UNUSED& ddx)
{return AUTO_HESS_EXT_VEC<TV,VEC,DIFF_UNUSED,1>(x,dx,ddx);}

template<class TV> AUTO_HESS_EXT_VEC<TV,DIFF_UNUSED,DIFF_UNUSED,0>
Make_Hess_Vec(const TV& x,const DIFF_UNUSED& dx,const DIFF_UNUSED& ddx)
{return AUTO_HESS_EXT_VEC<TV,DIFF_UNUSED,DIFF_UNUSED,0>(x,dx,ddx);}

template<class LAYOUT,class T,int i>
struct HESS_FROM_VAR_HELPER
{
    static_assert(is_scalar<T>::value,"This version of HESS_FROM_VAR_HELPER is only for scalars");

    typedef AUTO_HESS_EXT<T,typename ONE_NONZERO_VECTOR<1,i,LAYOUT,-1>::TYPE,typename EMPTY_MAT<LAYOUT,-1>::TYPE,3> TYPE;
};

template<class LAYOUT,class T,int d,int i>
struct HESS_FROM_VAR_HELPER<LAYOUT,VECTOR<T,d>,i>
{
    typedef AUTO_HESS_EXT_VEC<VECTOR<T,d>,typename ONE_NONZERO_VECTOR<1,i,LAYOUT,d>::TYPE,typename EMPTY_MAT<LAYOUT,d>::TYPE,3> TYPE;
};

template<class LAYOUT,int i,class T> inline
typename HESS_FROM_VAR_HELPER<LAYOUT,T,i>::TYPE Hess_From_Var(const T& v)
{return typename HESS_FROM_VAR_HELPER<LAYOUT,T,i>::TYPE(v);}

template<class LAYOUT,class T,int i>
struct DIFF_FROM_VAR_HELPER
{
    static_assert(is_scalar<T>::value,"This version of DIFF_FROM_VAR_HELPER is only for scalars");

    typedef AUTO_HESS_EXT<T,typename ONE_NONZERO_VECTOR<1,i,LAYOUT,-1>::TYPE,DIFF_UNUSED,1> TYPE;
};

template<class LAYOUT,class T,int d,int i>
struct DIFF_FROM_VAR_HELPER<LAYOUT,VECTOR<T,d>,i>
{
    typedef AUTO_HESS_EXT_VEC<VECTOR<T,d>,typename ONE_NONZERO_VECTOR<1,i,LAYOUT,d>::TYPE,DIFF_UNUSED,1> TYPE;
};

template<class LAYOUT,int i,class T> inline
typename DIFF_FROM_VAR_HELPER<LAYOUT,T,i>::TYPE Diff_From_Var(const T& v)
{return typename DIFF_FROM_VAR_HELPER<LAYOUT,T,i>::TYPE(v);}

template<class LAYOUT,int i,class T> inline
AUTO_HESS_EXT<T,DIFF_UNUSED,DIFF_UNUSED,0> From_Var(const T& v)
{return AUTO_HESS_EXT<T,DIFF_UNUSED,DIFF_UNUSED,0>(v);}

template<class LAYOUT,int i,class T,int d> inline
AUTO_HESS_EXT_VEC<VECTOR<T,d>,DIFF_UNUSED,DIFF_UNUSED,0> From_Var(const VECTOR<T,d>& v)
{return AUTO_HESS_EXT_VEC<VECTOR<T,d>,DIFF_UNUSED,DIFF_UNUSED,0>(v);}

template<class LAYOUT,class T,int d>
AUTO_HESS_EXT_VEC<VECTOR<T,d>,typename EMPTY_VEC<1,LAYOUT,d>::TYPE,typename EMPTY_MAT<LAYOUT,d>::TYPE,3>
Hess_From_Const(VECTOR<T,d> a)
{return AUTO_HESS_EXT_VEC<VECTOR<T,d>,typename EMPTY_VEC<1,LAYOUT,d>::TYPE,typename EMPTY_MAT<LAYOUT,d>::TYPE,3>(a);}

template<class LAYOUT,class T,int d>
AUTO_HESS_EXT_VEC<VECTOR<T,d>,typename EMPTY_VEC<1,LAYOUT,d>::TYPE,DIFF_UNUSED,1>
Diff_From_Const(VECTOR<T,d> a)
{return AUTO_HESS_EXT_VEC<VECTOR<T,d>,typename EMPTY_VEC<1,LAYOUT,d>::TYPE,DIFF_UNUSED,1>(a);}

template<class LAYOUT,class T,int d>
AUTO_HESS_EXT_VEC<VECTOR<T,d>,DIFF_UNUSED,DIFF_UNUSED,0>
From_Const(VECTOR<T,d> a)
{return AUTO_HESS_EXT_VEC<VECTOR<T,d>,DIFF_UNUSED,DIFF_UNUSED,0>(a);}

template<class T,int d,class VEC,class MAT,int Q> inline auto
operator*(const AUTO_HESS_EXT<T,VEC,MAT,Q>& a,const VECTOR<T,d>& v)
    -> decltype(Make_Hess_Vec(v*a.x,Outer_Product_Q<Q&1>(v,a.dx),Tensor_Product_0(a.ddx,v)))
{return Make_Hess_Vec(v*a.x,Outer_Product_Q<Q&1>(v,a.dx),Tensor_Product_0(a.ddx,v));}

template<class T,int d,class VEC,class MAT,int Q> inline auto
operator*(const VECTOR<T,d>& v,const AUTO_HESS_EXT<T,VEC,MAT,Q>& a) -> decltype(a*v)
{return a*v;}

template<class T,class TV,class VEC,class MAT,class VEC1,class MAT1,int Q> inline auto
operator*(const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& a,const AUTO_HESS_EXT_VEC<TV,VEC,MAT,Q>& v) -> decltype(v*a)
{return v*a;}

template<class TV,class VEC,class MAT,int Q> auto
operator*(typename TV::SCALAR a,const AUTO_HESS_EXT_VEC<TV,VEC,MAT,Q>& v) -> decltype(v*a)
{return v*a;}

template<class TV,class VEC,class MAT,int Q> auto
operator+(TV a,const AUTO_HESS_EXT_VEC<TV,VEC,MAT,Q>& u) -> decltype(Make_Hess_Vec(u.x+a,u.dx,u.ddx))
{return Make_Hess_Vec(u.x+a,u.dx,u.ddx);}

template<class TV,class VEC,class MAT,int Q> auto
operator-(TV a,const AUTO_HESS_EXT_VEC<TV,VEC,MAT,Q>& u)
    -> decltype(Make_Hess_Vec(a-u.x,-u.dx,-u.ddx))
{return Make_Hess_Vec(a-u.x,-u.dx,-u.ddx);}

template<class T,int d,class VEC1,class MAT1,int Q> auto
operator/(const VECTOR<T,d>& u,const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& a) -> decltype(u*((T)1/a))
{return u*((T)1/a);}

template<class TV,class VEC,class MAT,int Q> auto
operator*(const MATRIX<typename TV::SCALAR,TV::m>& m,const AUTO_HESS_EXT_VEC<TV,VEC,MAT,Q>& v)
    -> decltype(Make_Hess_Vec(m*v.x,m*v.dx,Contract_00(v.ddx,m.Transposed())))
{return Make_Hess_Vec(m*v.x,m*v.dx,Contract_00(v.ddx,m.Transposed()));}

template<class TV,class VEC,class MAT,int Q> inline auto
sin(const AUTO_HESS_EXT_VEC<TV,VEC,MAT,Q>& a)
    -> decltype(Make_Hess_Vec(TV(),DIAGONAL_MATRIX<typename TV::SCALAR,TV::m>()*a.dx,Contract_00(a.ddx,DIAGONAL_MATRIX<typename TV::SCALAR,TV::m>())-Symmetric_Double_Contract_12_With_Tensor_Q<Q&2>(DIAGONAL_TENSOR<typename TV::SCALAR,TV::m>(),a.dx,a.dx)))
{TV s=sin(a.x),c=cos(a.x);return Make_Hess_Vec(s,DIAGONAL_MATRIX<typename TV::SCALAR,TV::m>(c)*a.dx,Contract_00(a.ddx,DIAGONAL_MATRIX<typename TV::SCALAR,TV::m>(c))-Symmetric_Double_Contract_12_With_Tensor_Q<Q&2>(DIAGONAL_TENSOR<typename TV::SCALAR,TV::m>(s/2),a.dx,a.dx));}

template<class TV,class VEC,class MAT,int Q> inline auto
cos(const AUTO_HESS_EXT_VEC<TV,VEC,MAT,Q>& a)
    -> decltype(Make_Hess_Vec(TV(),-DIAGONAL_MATRIX<typename TV::SCALAR,TV::m>()*a.dx,Contract_00(a.ddx,-DIAGONAL_MATRIX<typename TV::SCALAR,TV::m>())-Symmetric_Double_Contract_12_With_Tensor_Q<Q&2>(DIAGONAL_TENSOR<typename TV::SCALAR,TV::m>(),a.dx,a.dx)))
{TV s=sin(a.x),c=cos(a.x);return Make_Hess_Vec(c,-DIAGONAL_MATRIX<typename TV::SCALAR,TV::m>(s)*a.dx,Contract_00(a.ddx,-DIAGONAL_MATRIX<typename TV::SCALAR,TV::m>(s))-Symmetric_Double_Contract_12_With_Tensor_Q<Q&2>(DIAGONAL_TENSOR<typename TV::SCALAR,TV::m>(c/2),a.dx,a.dx));}

template<class T_MAT,class VEC,class MAT,int Q> struct AUTO_HESS_EXT_MAT;

template<class T_MAT,class VEC,class MAT> AUTO_HESS_EXT_MAT<T_MAT,VEC,MAT,3>
Make_Hess_Mat(const T_MAT& x,const GRADIENT_MAT<MATRIX<typename T_MAT::SCALAR,T_MAT::m,T_MAT::n>,VEC>& dx,const HESSIAN_MAT<MATRIX<typename T_MAT::SCALAR,T_MAT::m,T_MAT::n>,MAT>& ddx);

template<class T_MAT,class VEC> AUTO_HESS_EXT_MAT<T_MAT,VEC,DIFF_UNUSED,1>
Make_Hess_Mat(const T_MAT& x,const GRADIENT_MAT<MATRIX<typename T_MAT::SCALAR,T_MAT::m,T_MAT::n>,VEC>& dx,const DIFF_UNUSED& ddx);

template<class T_MAT> AUTO_HESS_EXT_MAT<T_MAT,DIFF_UNUSED,DIFF_UNUSED,0>
Make_Hess_Mat(const T_MAT& x,const DIFF_UNUSED& dx,const DIFF_UNUSED& ddx);

template<class T_MAT,class VEC,class MAT,int Q>
struct AUTO_HESS_EXT_MAT
{
    typedef typename T_MAT::SCALAR T;
    typedef typename conditional<Q&1,GRADIENT_MAT<MATRIX<typename T_MAT::SCALAR,T_MAT::m,T_MAT::n>,VEC>,DIFF_UNUSED>::type GRAD;
    typedef typename conditional<Q&2,HESSIAN_MAT<MATRIX<typename T_MAT::SCALAR,T_MAT::m,T_MAT::n>,MAT>,DIFF_UNUSED>::type HESS;

    T_MAT x;
    GRAD dx;
    HESS ddx;

    AUTO_HESS_EXT_MAT(T_MAT x=T_MAT()):
        x(x)
    {}

    AUTO_HESS_EXT_MAT(T_MAT x,const GRAD& dx,const HESS& ddx):
        x(x),dx(dx),ddx(ddx)
    {}

    auto operator-() const -> decltype(Make_Hess_Mat(-this->x,-this->dx,-this->ddx))
    {return Make_Hess_Mat(-x,-dx,-ddx);}

    AUTO_HESS_EXT_MAT operator+() const
    {return *this;}

    template<class T_MAT1,class VEC1,class MAT1> auto
    operator+(const AUTO_HESS_EXT_MAT<T_MAT1,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess_Mat(x+a.x,dx+a.dx,ddx+a.ddx))
    {return Make_Hess_Mat(x+a.x,dx+a.dx,ddx+a.ddx);}

    template<class T_MAT1,class VEC1,class MAT1> auto
    operator-(const AUTO_HESS_EXT_MAT<T_MAT1,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess_Mat(x-a.x,dx-a.dx,ddx-a.ddx))
    {return Make_Hess_Mat(x-a.x,dx-a.dx,ddx-a.ddx);}

    template<class VEC1,class MAT1> auto
    operator*(const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess_Mat(x*a.x,a.x*dx+Tensor_Product_2_Q<Q&1>(x,a.dx),ddx*a.x+Tensor_Product_01(x,a.ddx)+Symmetric_Tensor_Product_23_Q<Q&2>(dx,a.dx)))
    {return Make_Hess_Mat(x*a.x,a.x*dx+Tensor_Product_2_Q<Q&1>(x,a.dx),ddx*a.x+Tensor_Product_01(x,a.ddx)+Symmetric_Tensor_Product_23_Q<Q&2>(dx,a.dx));}

    template<class VEC1,class MAT1> auto
    operator/(const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& a) const
        -> decltype(*this*((T)1/a))
    {return *this*((T)1/a);}

    template<class T_MAT2>
    auto operator+(const T_MAT2& a) const
        -> typename enable_if<IS_MATRIX<T_MAT2>::value,decltype(Make_Hess_Mat(this->x+a,this->dx,this->ddx))>::type
    {return Make_Hess_Mat(x+a,dx,ddx);}

    template<class T_MAT2>
    auto operator-(const T_MAT2& a) const
        -> typename enable_if<IS_MATRIX<T_MAT2>::value,decltype(Make_Hess_Mat(this->x-a,this->dx,this->ddx))>::type
    {return Make_Hess_Mat(x-a,dx,ddx);}

    auto operator*(T a) const -> decltype(Make_Hess_Mat(this->x*a,this->dx*a,this->ddx*a))
    {return Make_Hess_Mat(x*a,dx*a,ddx*a);}

    auto operator/(T a) const -> decltype(Make_Hess_Mat(this->x/a,this->dx/a,this->ddx/a))
    {return Make_Hess_Mat(x/a,dx/a,ddx/a);}
    
    template<class TV,class VEC1,class MAT1>
    auto Transpose_Times(AUTO_HESS_EXT_VEC<TV,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess_Vec(Transpose_Times(this->x,a.x),Contract_0(this->dx,a.x)+Transpose_Times(this->x,a.dx),Contract_0(this->ddx,a.x)+Symmetric_Contract_00_12_Q<Q&2>(this->dx,a.dx)+Contract_00(a.ddx,this->x)))
    {
        return Make_Hess_Vec(Transpose_Times(x,a.x),Contract_0(dx,a.x)+Transpose_Times(x,a.dx),Contract_0(ddx,a.x)+Symmetric_Contract_00_12_Q<Q&2>(dx,a.dx)+Contract_00(a.ddx,x));
    }

    template<int d>
    auto Transpose_Times(const VECTOR<T,d>& a) const
        -> decltype(Make_Hess_Vec(Transpose_Times(this->x,a),Contract_0(this->dx,a),Contract_0(this->ddx,a)))
    {
        return Make_Hess_Vec(Transpose_Times(x,a),Contract_0(dx,a),Contract_0(ddx,a));
    }

    auto Trace() const -> decltype(Make_Hess(this->x.Trace(),Contract_01_Q<Q&1>(this->dx),Contract_01_Q<Q&2>(this->ddx)))
    {return Make_Hess(x.Trace(),Contract_01_Q<Q&1>(dx),Contract_01_Q<Q&2>(ddx));}

    auto Transposed() const -> decltype(Make_Hess_Mat(this->x.Transposed(),Transposed_01_Q<Q&1>(this->dx),Transposed_01_Q<Q&2>(this->ddx)))
    {return Make_Hess_Mat(x.Transposed(),Transposed_01_Q<Q&1>(dx),Transposed_01_Q<Q&2>(ddx));}

    auto Twice_Symmetric_Part() const -> decltype(Make_Hess_Mat(this->x.Twice_Symmetric_Part(),Twice_Symmetric_Part_01_Q<Q&1>(this->dx),Twice_Symmetric_Part_01_Q<Q&2>(this->ddx)))
    {return Make_Hess_Mat(x.Twice_Symmetric_Part(),Twice_Symmetric_Part_01_Q<Q&1>(dx),Twice_Symmetric_Part_01_Q<Q&2>(ddx));}

    template<class T_MAT1,class VEC1,class MAT1> auto
    operator*(const AUTO_HESS_EXT_MAT<T_MAT1,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess_Mat(x*a.x,Contract_10(dx,a.x)+Contract_01_Q<Q&1>(a.dx,x),Contract_10(ddx,a.x)+Contract_01_Q<Q&2>(a.ddx,x)+Symmetric_Contract_01_23_Q<Q&2>(a.dx,dx)))
    {return Make_Hess_Mat(x*a.x,Contract_10(dx,a.x)+Contract_01_Q<Q&1>(a.dx,x),Contract_10(ddx,a.x)+Contract_01_Q<Q&2>(a.ddx,x)+Symmetric_Contract_01_23_Q<Q&2>(a.dx,dx));}

    template<class T_MAT1,class VEC1,class MAT1> auto
    Transpose_Times(const AUTO_HESS_EXT_MAT<T_MAT1,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess_Mat(Transpose_Times(x,a.x),Contract_00(dx,a.x)+Contract_00(a.dx,x),Contract_00(ddx,a.x)+Contract_00(a.ddx,x)+Symmetric_Contract_00_23_Q<Q&2>(dx,a.dx)))
    {return Make_Hess_Mat(Transpose_Times(x,a.x),Contract_00(dx,a.x)+Contract_00(a.dx,x),Contract_00(ddx,a.x)+Contract_00(a.ddx,x)+Symmetric_Contract_00_23_Q<Q&2>(dx,a.dx));}

    template<class T_MAT1,class VEC1,class MAT1> auto
    Times_Transpose(const AUTO_HESS_EXT_MAT<T_MAT1,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess_Mat(x*a.x,Contract_11(dx,a.x)+Contract_11(a.dx,x),Contract_11(ddx,a.x)+Contract_11(a.ddx,x)+Symmetric_Contract_11_23_Q<Q&2>(dx,a.dx)))
    {return Make_Hess_Mat(x*a.x,Contract_11(dx,a.x)+Contract_11(a.dx,x),Contract_11(ddx,a.x)+Contract_11(a.ddx,x)+Symmetric_Contract_11_23_Q<Q&2>(dx,a.dx));}

    template<class T_MAT1> auto
    operator*(const T_MAT1& a) const
        -> decltype(Make_Hess_Mat(x*a,Contract_10(dx,a),Contract_10(ddx,a)))
    {return Make_Hess_Mat(x*a,Contract_10(dx,a),Contract_10(ddx,a));}

    template<class T_MAT1> auto
    Transpose_Times(const T_MAT1& a) const
        -> decltype(Make_Hess_Mat(Transpose_Times(x,a),Contract_00(dx,a),Contract_00(ddx,a)))
    {return Make_Hess_Mat(Transpose_Times(x,a),Contract_00(dx,a),Contract_00(ddx,a));}

    template<class T_MAT1> auto
    Times_Transpose(const T_MAT1& a) const
        -> decltype(Make_Hess_Mat(Times_Transpose(x,a),Contract_11(dx,a),Contract_11(ddx,a)))
    {return Make_Hess_Mat(Times_Transpose(x,a),Contract_11(dx,a),Contract_11(ddx,a));}

    auto Normal_Equations_Matrix() const
        -> decltype(Make_Hess_Mat(x.Normal_Equations_Matrix(),Twice_Symmetric_Part_01_Q<Q&1>(Contract_00(dx,x)),Twice_Symmetric_Part_01_Q<Q&2>(Contract_00(ddx,x))+Symmetric_Contract_00_23_Q<Q&2>(dx,dx)))
    {return Make_Hess_Mat(x.Normal_Equations_Matrix(),Twice_Symmetric_Part_01_Q<Q&1>(Contract_00(dx,x)),Twice_Symmetric_Part_01_Q<Q&2>(Contract_00(ddx,x))+Symmetric_Contract_00_23_Q<Q&2>(dx,dx));}

    auto Outer_Product_Matrix() const
        -> decltype(Make_Hess_Mat(x.Outer_Product_Matrix(),Twice_Symmetric_Part_01_Q<Q&1>(Contract_11(dx,x)),Twice_Symmetric_Part_01_Q<Q&2>(Contract_11(ddx,x))+Symmetric_Contract_11_23_Q<Q&2>(dx,dx)))
    {return Make_Hess_Mat(x.Outer_Product_Matrix(),Twice_Symmetric_Part_01_Q<Q&1>(Contract_11(dx,x)),Twice_Symmetric_Part_01_Q<Q&2>(Contract_11(ddx,x))+Symmetric_Contract_11_23_Q<Q&2>(dx,dx));}

    template<class T_MAT1,class VEC1,class MAT1> auto
    Double_Contract(const AUTO_HESS_EXT_MAT<T_MAT1,VEC1,MAT1,Q>& a) const
        -> decltype(Make_Hess(Double_Contract(x,a.x),Double_Contract_00_11(dx,a.x)+Double_Contract_00_11(a.dx,x),Double_Contract_00_11(ddx,a.x)+Double_Contract_00_11(a.ddx,x)+Symmetric_Double_Contract_00_11_01_Q<Q&2>(a.dx,dx)))
    {return Make_Hess(Double_Contract(x,a.x),Double_Contract_00_11(dx,a.x)+Double_Contract_00_11(a.dx,x),Double_Contract_00_11(ddx,a.x)+Double_Contract_00_11(a.ddx,x)+Symmetric_Double_Contract_00_11_01_Q<Q&2>(a.dx,dx));}

    template<class T_MAT1> auto
    Double_Contract(const T_MAT1& a) const
        -> decltype(Make_Hess(Double_Contract(a,x),Double_Contract_00_11(dx,a),Double_Contract_00_11(ddx,a)))
    {return Make_Hess(Double_Contract(x,a),Double_Contract_00_11(dx,a),Double_Contract_00_11(ddx,a));}
};

template<class T_MAT,class VEC,class MAT> AUTO_HESS_EXT_MAT<T_MAT,VEC,MAT,3>
Make_Hess_Mat(const T_MAT& x,const GRADIENT_MAT<MATRIX<typename T_MAT::SCALAR,T_MAT::m,T_MAT::n>,VEC>& dx,const HESSIAN_MAT<MATRIX<typename T_MAT::SCALAR,T_MAT::m,T_MAT::n>,MAT>& ddx)
{return AUTO_HESS_EXT_MAT<T_MAT,VEC,MAT,3>(x,dx,ddx);}

template<class T_MAT,class VEC> AUTO_HESS_EXT_MAT<T_MAT,VEC,DIFF_UNUSED,1>
Make_Hess_Mat(const T_MAT& x,const GRADIENT_MAT<MATRIX<typename T_MAT::SCALAR,T_MAT::m,T_MAT::n>,VEC>& dx,const DIFF_UNUSED& ddx)
{return AUTO_HESS_EXT_MAT<T_MAT,VEC,DIFF_UNUSED,1>(x,dx,ddx);}

template<class T_MAT> AUTO_HESS_EXT_MAT<T_MAT,DIFF_UNUSED,DIFF_UNUSED,0>
Make_Hess_Mat(const T_MAT& x,const DIFF_UNUSED& dx,const DIFF_UNUSED& ddx)
{return AUTO_HESS_EXT_MAT<T_MAT,DIFF_UNUSED,DIFF_UNUSED,0>(x,dx,ddx);}

template<class LAYOUT,class A>
typename enable_if<IS_MATRIX<A>::value,AUTO_HESS_EXT_MAT<A,typename EMPTY_VEC<1,LAYOUT,A::m,A::n>::TYPE,typename EMPTY_MAT<LAYOUT,A::m,A::n>::TYPE,3> >::type
Hess_From_Const(const A& a)
{return AUTO_HESS_EXT_MAT<A,typename EMPTY_VEC<1,LAYOUT,A::m,A::n>::TYPE,typename EMPTY_MAT<LAYOUT,A::m,A::n>::TYPE,3>(a);}

template<class LAYOUT,class A>
typename enable_if<IS_MATRIX<A>::value,AUTO_HESS_EXT_MAT<A,typename EMPTY_VEC<1,LAYOUT,A::m,A::n>::TYPE,DIFF_UNUSED,1> >::type
Diff_From_Const(const A& a)
{return AUTO_HESS_EXT_MAT<A,typename EMPTY_VEC<1,LAYOUT,A::m,A::n>::TYPE,DIFF_UNUSED,1>(a);}

template<class LAYOUT,class A>
typename enable_if<IS_MATRIX<A>::value,AUTO_HESS_EXT_MAT<A,DIFF_UNUSED,DIFF_UNUSED,0> >::type
From_Const(const A& a)
{return AUTO_HESS_EXT_MAT<A,DIFF_UNUSED,DIFF_UNUSED,0>(a);}

template<class T_MAT2,class T_MAT,class VEC,class MAT,int Q>
auto operator+(const T_MAT2& a,const AUTO_HESS_EXT_MAT<T_MAT,VEC,MAT,Q>& b)
    -> typename enable_if<IS_MATRIX<T_MAT2>::value,decltype(Make_Hess_Mat(a+b.x,b.dx,b.ddx))>::type
{return Make_Hess_Mat(a+b.x,b.dx,b.ddx);}

template<class T_MAT2,class T_MAT,class VEC,class MAT,int Q>
auto operator-(const T_MAT2& a,const AUTO_HESS_EXT_MAT<T_MAT,VEC,MAT,Q>& b)
    -> typename enable_if<IS_MATRIX<T_MAT2>::value,decltype(Make_Hess_Mat(a-b.x,-b.dx,-b.ddx))>::type
{return Make_Hess_Mat(a-b.x,-b.dx,-b.ddx);}

template<class T_MAT,class VEC,class MAT,int Q>
auto operator*(const typename T_MAT::SCALAR& a,const AUTO_HESS_EXT_MAT<T_MAT,VEC,MAT,Q>& b) -> decltype(b*a)
{return b*a;}

template<class T,class VEC1,class MAT1,class T_MAT,class VEC,class MAT,int Q>
auto operator*(const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& a,const AUTO_HESS_EXT_MAT<T_MAT,VEC,MAT,Q>& b) -> decltype(b*a)
{return b*a;}

template<class TV,class VEC1,class MAT1,class T_MAT,class VEC,class MAT,int Q>
auto operator*(const AUTO_HESS_EXT_MAT<T_MAT,VEC,MAT,Q>& a,const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1,Q>& b)
    -> decltype(Make_Hess_Vec(a.x*b.x,Contract_1(a.dx,b.x)+a.x*b.dx,Contract_1(a.ddx,b.x)+Symmetric_Contract_10_12_Q<Q&2>(a.dx,b.dx)+Contract_01_Q<Q&2>(b.ddx,a.x)))
{
    return Make_Hess_Vec(a.x*b.x,Contract_1(a.dx,b.x)+a.x*b.dx,Contract_1(a.ddx,b.x)+Symmetric_Contract_10_12_Q<Q&2>(a.dx,b.dx)+Contract_01_Q<Q&2>(b.ddx,a.x));
}

template<class T,int d,class T_MAT,class VEC,class MAT,int Q>
auto operator*(const AUTO_HESS_EXT_MAT<T_MAT,VEC,MAT,Q>& a,const VECTOR<T,d>& b)
    -> decltype(Make_Hess_Vec(a.x*b,Contract_1(a.dx,b),Contract_1(a.ddx,b)))
{
    return Make_Hess_Vec(a.x*b,Contract_1(a.dx,b),Contract_1(a.ddx,b));
}

template<class TV,class VEC1,class MAT1,class T_MAT,int Q>
auto operator*(const T_MAT& a,const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1,Q>& b)
    -> typename enable_if<IS_MATRIX<T_MAT>::value,decltype(Make_Hess_Vec(a*b.x,a*b.dx,Contract_01_Q<Q&2>(b.ddx,a)))>::type
{
    return Make_Hess_Vec(a*b.x,a*b.dx,Contract_01_Q<Q&2>(b.ddx,a));
}

template<class T,class VEC1,class MAT1,class T_MAT,int Q>
auto operator*(const T_MAT& a,const AUTO_HESS_EXT<T,VEC1,MAT1,Q>& b)
    -> typename enable_if<IS_MATRIX<T_MAT>::value,decltype(Make_Hess_Mat(a*b.x,Tensor_Product_2_Q<Q&1>(a,b.dx),Tensor_Product_01(a,b.ddx)))>::type
{
    return Make_Hess_Mat(a*b.x,Tensor_Product_2_Q<Q&1>(a,b.dx),Tensor_Product_01(a,b.ddx));
}

template<class T_MAT,class VEC,class MAT,class TV,int Q>
auto Transpose_Times(const AUTO_HESS_EXT_MAT<T_MAT,VEC,MAT,Q>& a,TV& b)
    -> decltype(Make_Hess_Vec(a.x.Transpose_Times(b),Contract_0(a.dx,b),Contract_0(a.ddx,b)))
{
    return Make_Hess_Vec(a.x.Transpose_Times(b),Contract_0(a.dx,b),Contract_0(a.ddx,b));
}

template<class T_MAT,class VEC,class MAT,class T_MAT1,int Q> auto
Transpose_Times(const AUTO_HESS_EXT_MAT<T_MAT,VEC,MAT,Q>& a,const T_MAT1& b)
    -> decltype(Make_Hess_Mat(a.x.Transpose_Times(b),Contract_00(a.dx,b),Contract_00(a.ddx,b)))
{return Make_Hess_Mat(a.x.Transpose_Times(b),Contract_00(a.dx,b),Contract_00(a.ddx,b));}

template<class T_MAT,class VEC,class MAT,class T_MAT1,int Q> auto
Times_Transpose(const AUTO_HESS_EXT_MAT<T_MAT,VEC,MAT,Q>& a,const T_MAT1& b)
    -> decltype(Make_Hess_Mat(a.x*b,Contract_11(a.dx,b),Contract_11(a.ddx,b)))
{return Make_Hess_Mat(a.x.Times_Transpose(b),Contract_11(a.dx,b),Contract_11(a.ddx,b));}

template<class TV,class VEC,class MAT,class VEC1,class MAT1,int Q> auto
Outer_Product(const AUTO_HESS_EXT_VEC<TV,VEC,MAT,Q>& a,const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1,Q>& b)
    -> decltype(Make_Hess_Mat(Outer_Product(a.x,b.x),Tensor_Product_1_Q<Q&1>(a.dx,b.x)+Tensor_Product_0(b.dx,a.x),Tensor_Product_1_Q<Q&2>(a.ddx,b.x)+Tensor_Product_0(b.ddx,a.x)+Symmetric_Tensor_Product_02_23_Q<Q&2>(a.dx,b.dx)))
{return Make_Hess_Mat(Outer_Product(a.x,b.x),Tensor_Product_1_Q<Q&1>(a.dx,b.x)+Tensor_Product_0(b.dx,a.x),Tensor_Product_1_Q<Q&2>(a.ddx,b.x)+Tensor_Product_0(b.ddx,a.x)+Symmetric_Tensor_Product_02_23_Q<Q&2>(a.dx,b.dx));}

template<class TV,class VEC1,class MAT1,int Q> auto
Outer_Product(const TV& a,const AUTO_HESS_EXT_VEC<TV,VEC1,MAT1,Q>& b)
    -> decltype(Make_Hess_Mat(Outer_Product(a,b.x),Tensor_Product_0(b.dx,a),Tensor_Product_0(b.ddx,a)))
{return Make_Hess_Mat(Outer_Product(a,b.x),Tensor_Product_0(b.dx,a),Tensor_Product_0(b.ddx,a));}

template<class TV,class VEC,class MAT,int Q> auto
Outer_Product(const AUTO_HESS_EXT_VEC<TV,VEC,MAT,Q>& a,const TV& b)
    -> decltype(Make_Hess_Mat(Outer_Product(a.x,b),Tensor_Product_1_Q<Q&1>(a.dx,b),Tensor_Product_1_Q<Q&2>(a.ddx,b)))
{return Make_Hess_Mat(Outer_Product(a.x,b),Tensor_Product_1_Q<Q&1>(a.dx,b),Tensor_Product_1_Q<Q&2>(a.ddx,b));}

template<class TV,class VEC,class MAT,int Q> auto
Outer_Product(const AUTO_HESS_EXT_VEC<TV,VEC,MAT,Q>& a)
    -> decltype(Make_Hess_Mat(Outer_Product(a.x),Twice_Symmetric_Part_01_Q<Q&1>(Tensor_Product_1_Q<Q&1>(a.dx,a.x)),Twice_Symmetric_Part_01_Q<Q&2>(Tensor_Product_1_Q<Q&2>(a.ddx,a.x))+Symmetric_Tensor_Product_02_23_Q<Q&2>(a.dx,a.dx)))
{return Make_Hess_Mat(Outer_Product(a.x),Twice_Symmetric_Part_01_Q<Q&1>(Tensor_Product_1_Q<Q&1>(a.dx,a.x)),Twice_Symmetric_Part_01_Q<Q&2>(Tensor_Product_1_Q<Q&2>(a.ddx,a.x))+Symmetric_Tensor_Product_02_23_Q<Q&2>(a.dx,a.dx));}
}
using HETERO_DIFF::AUTO_HESS_EXT;
using HETERO_DIFF::AUTO_HESS_EXT_VEC;
using HETERO_DIFF::AUTO_HESS_EXT_MAT;
using HETERO_DIFF::Hess_From_Const;
using HETERO_DIFF::Hess_From_Var;
using HETERO_DIFF::Diff_From_Const;
using HETERO_DIFF::Diff_From_Var;
using HETERO_DIFF::From_Const;
using HETERO_DIFF::From_Var;
using HETERO_DIFF::Extract;
using HETERO_DIFF::Get;
using HETERO_DIFF::sin;
using HETERO_DIFF::cos;
using HETERO_DIFF::tan;
//using HETERO_DIFF::atan2;
using HETERO_DIFF::exp;
using HETERO_DIFF::log;
// using HETERO_DIFF::min;
// using HETERO_DIFF::max;
using HETERO_DIFF::abs;
//using HETERO_DIFF::hypot;
using HETERO_DIFF::sqr;
using HETERO_DIFF::sqrt;
}
#endif
