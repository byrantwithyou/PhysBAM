//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_SIGNED_VOLUME__
#define __ADAPTIVE_SIGNED_VOLUME__
#include <Tools/Adaptive_Arithmetic/ADAPTIVE_DETERMINANT.h>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{

template<class T_ADAPTIVE,int d> struct IS_ADAPTIVE<ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,d> > {static const bool value=true;};

template<class T_ADAPTIVE>
class ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,1>
    :public ADAPTIVE_WRAPPER_BASE<ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,1>,ADAPTIVE_DETERMINANT<ADAPTIVE_DIFFERENCE<T_ADAPTIVE,T_ADAPTIVE>,1> >
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
public:
    typedef typename ADAPTIVE_WRAPPER_BASE<ADAPTIVE_SIGNED_VOLUME,ADAPTIVE_DETERMINANT<ADAPTIVE_DIFFERENCE<T_ADAPTIVE,T_ADAPTIVE>,1> >::WRAPPED_TYPE WRAPPED_TYPE;
private:
    const T_ADAPTIVE a00,a10;
    ADAPTIVE_SIGNED_VOLUME();
    ADAPTIVE_SIGNED_VOLUME& operator=(const ADAPTIVE_SIGNED_VOLUME&);

public:
    ADAPTIVE_SIGNED_VOLUME(const T_ADAPTIVE& a11_input,const T_ADAPTIVE& a21_input)
        :a00(a11_input),a10(a21_input)
    {}

public:
    WRAPPED_TYPE Wrapped_Implementation() const
    {return Adaptive_Determinant<void>(a00-a10);}
//#####################################################################
};

template<class T_ADAPTIVE>
class ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,2>
    :public ADAPTIVE_WRAPPER_BASE<ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,2>,ADAPTIVE_DETERMINANT<ADAPTIVE_DIFFERENCE<T_ADAPTIVE,T_ADAPTIVE>,2> >
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
public:
    typedef typename ADAPTIVE_WRAPPER_BASE<ADAPTIVE_SIGNED_VOLUME,ADAPTIVE_DETERMINANT<ADAPTIVE_DIFFERENCE<T_ADAPTIVE,T_ADAPTIVE>,2> >::WRAPPED_TYPE WRAPPED_TYPE;
private:
    const T_ADAPTIVE a00,a01,a10,a11,a20,a21;
    ADAPTIVE_SIGNED_VOLUME();
    ADAPTIVE_SIGNED_VOLUME& operator=(const ADAPTIVE_SIGNED_VOLUME&);

public:
    ADAPTIVE_SIGNED_VOLUME(const T_ADAPTIVE& a11_input,const T_ADAPTIVE& a12_input,const T_ADAPTIVE& a21_input,const T_ADAPTIVE& a22_input,const T_ADAPTIVE& a31_input,const T_ADAPTIVE& a32_input)
        :a00(a11_input),a01(a12_input),a10(a21_input),a11(a22_input),a20(a31_input),a21(a32_input)
    {}

public:
    WRAPPED_TYPE Wrapped_Implementation() const
    {return Adaptive_Determinant<void>(a00-a20,a01-a21,a10-a20,a11-a21);}
//#####################################################################
};

template<class T_ADAPTIVE>
class ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,3>
    :public ADAPTIVE_WRAPPER_BASE<ADAPTIVE_SIGNED_VOLUME<T_ADAPTIVE,3>,ADAPTIVE_DETERMINANT<ADAPTIVE_DIFFERENCE<T_ADAPTIVE,T_ADAPTIVE>,3> >
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
public:
    typedef typename ADAPTIVE_WRAPPER_BASE<ADAPTIVE_SIGNED_VOLUME,ADAPTIVE_DETERMINANT<ADAPTIVE_DIFFERENCE<T_ADAPTIVE,T_ADAPTIVE>,3> >::WRAPPED_TYPE WRAPPED_TYPE;
private:
    const T_ADAPTIVE a00,a01,a02,a10,a11,a12,a20,a21,a22,a30,a31,a32;
    ADAPTIVE_SIGNED_VOLUME();
    ADAPTIVE_SIGNED_VOLUME& operator=(const ADAPTIVE_SIGNED_VOLUME&);

public:
    ADAPTIVE_SIGNED_VOLUME(const T_ADAPTIVE& a11_input,const T_ADAPTIVE& a12_input,const T_ADAPTIVE& a13_input,const T_ADAPTIVE& a21_input,const T_ADAPTIVE& a22_input,const T_ADAPTIVE& a23_input,
        const T_ADAPTIVE& a31_input,const T_ADAPTIVE& a32_input,const T_ADAPTIVE& a33_input,const T_ADAPTIVE& a41_input,const T_ADAPTIVE& a42_input,const T_ADAPTIVE& a43_input)
        :a00(a11_input),a01(a12_input),a02(a13_input),a10(a21_input),a11(a22_input),a12(a23_input),a20(a31_input),a21(a32_input),a22(a33_input),a30(a41_input),a31(a42_input),a32(a43_input)
    {}

public:
    WRAPPED_TYPE Wrapped_Implementation() const
    {return Adaptive_Determinant<void>(a00-a30,a01-a31,a02-a32,a10-a30,a11-a31,a12-a32,a20-a30,a21-a31,a22-a32);}
//#####################################################################
};

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Signed_Volume(const T& t00,const T& t10)
{return ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>(t00,t10);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Signed_Volume(const VECTOR<T,1>& t1,const VECTOR<T,1>& t2)
{return Adaptive_Signed_Volume<T_EXACT>(t1[0],t2[0]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Signed_Volume(const VECTOR<T,1>& t1,const VECTOR<VECTOR<T,1>,1>& t2)
{return Adaptive_Signed_Volume<T_EXACT>(t1,t2[0]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Signed_Volume(const VECTOR<VECTOR<T,1>,2>& t)
{return Adaptive_Signed_Volume<T_EXACT,T>(t[0],t[0]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Signed_Volume(const T& t00,const T& t01,const T& t10,const T& t11,const T& t20,const T& t21)
{return ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>(t00,t01,t10,t11,t20,t21);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Signed_Volume(const VECTOR<T,2>& t1,const VECTOR<T,2>& t2,const VECTOR<T,2>& t3)
{return Adaptive_Signed_Volume<T_EXACT>(t1[0],t1[0],t2[0],t2[0],t3[0],t3[0]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Signed_Volume(const VECTOR<T,2>& t1,const VECTOR<VECTOR<T,2>,2>& t12)
{return Adaptive_Signed_Volume<T_EXACT>(t1,t12[0],t12[0]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Signed_Volume(const VECTOR<VECTOR<T,2>,3>& t)
{return Adaptive_Signed_Volume<T_EXACT>(t[0],t[0],t[1]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Signed_Volume(const T& t00,const T& t01,const T& t02,const T& t10,const T& t11,const T& t12,const T& t20,const T& t21,const T& t22,const T& t30,const T& t31,const T& t32)
{return ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>(t00,t01,t02,t10,t11,t12,t20,t21,t22,t30,t31,t32);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Signed_Volume(const VECTOR<T,3>& t1,const VECTOR<T,3>& t2,const VECTOR<T,3>& t3,const VECTOR<T,3>& t4)
{return Adaptive_Signed_Volume<T_EXACT>(t1[0],t1[0],t1[1],t2[0],t2[0],t2[1],t3[0],t3[0],t3[1],t4[0],t4[0],t4[1]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Signed_Volume(const VECTOR<T,3>& t1,const VECTOR<VECTOR<T,3>,3>& t234)
{return Adaptive_Signed_Volume<T_EXACT>(t1,t234[0],t234[0],t234[1]);}

template<class T_EXACT,class T>
ADAPTIVE_SIGNED_VOLUME<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Signed_Volume(const VECTOR<VECTOR<T,3>,4>& t)
{return Adaptive_Signed_Volume<T_EXACT>(t[0],t[0],t[1],t[2]);}

}
using ADAPTIVE_DETAIL::ADAPTIVE_SIGNED_VOLUME;
using ADAPTIVE_DETAIL::Adaptive_Signed_Volume;
}
#endif
