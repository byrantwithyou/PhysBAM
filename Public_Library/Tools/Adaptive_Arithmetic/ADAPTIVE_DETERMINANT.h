//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_DETERMINANT__
#define __ADAPTIVE_DETERMINANT__
#include <Tools/Adaptive_Arithmetic/ADAPTIVE_FORWARD.h>
#include <Tools/Adaptive_Arithmetic/ADAPTIVE_OP.h>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{

template<class T_ADAPTIVE,int d> struct IS_ADAPTIVE<ADAPTIVE_DETERMINANT<T_ADAPTIVE,d> > {static const bool value=true;};

template<class T_ADAPTIVE>
class ADAPTIVE_DETERMINANT<T_ADAPTIVE,1>
    :public ADAPTIVE_WRAPPER_BASE<ADAPTIVE_DETERMINANT<T_ADAPTIVE,1>,T_ADAPTIVE>
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
public:
    typedef typename ADAPTIVE_WRAPPER_BASE<ADAPTIVE_DETERMINANT,T_ADAPTIVE>::WRAPPED_TYPE WRAPPED_TYPE;
private:
    const T_ADAPTIVE a00;
    ADAPTIVE_DETERMINANT();
    void operator=(const ADAPTIVE_DETERMINANT&);

public:
    ADAPTIVE_DETERMINANT(const T_ADAPTIVE& a11_input)
        :a00(a11_input)
    {}

    WRAPPED_TYPE Wrapped_Implementation() const
    {return a00;}
//#####################################################################
};

template<class T_ADAPTIVE>
class ADAPTIVE_DETERMINANT<T_ADAPTIVE,2>
    :public ADAPTIVE_WRAPPER_BASE<ADAPTIVE_DETERMINANT<T_ADAPTIVE,2>,typename ADAPTIVE_VECTOR_OPERATION_POLICY<NIL,T_ADAPTIVE,2>::DETERMINANT>
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
public:
    typedef typename ADAPTIVE_WRAPPER_BASE<ADAPTIVE_DETERMINANT,typename ADAPTIVE_VECTOR_OPERATION_POLICY<NIL,T_ADAPTIVE,2>::DETERMINANT>::WRAPPED_TYPE WRAPPED_TYPE;
private:
    const T_ADAPTIVE a00,a01,a10,a11;
    ADAPTIVE_DETERMINANT();
    void operator=(const ADAPTIVE_DETERMINANT&);

public:
    ADAPTIVE_DETERMINANT(const T_ADAPTIVE& a11_input,const T_ADAPTIVE& a12_input,const T_ADAPTIVE& a21_input,const T_ADAPTIVE& a22_input)
        :a00(a11_input),a01(a12_input),a10(a21_input),a11(a22_input)
    {}

    WRAPPED_TYPE Wrapped_Implementation() const
    {return (a00*a11-a01*a10);}
//#####################################################################
};

template<class T_ADAPTIVE>
class ADAPTIVE_DETERMINANT<T_ADAPTIVE,3>
    :public ADAPTIVE_WRAPPER_BASE<ADAPTIVE_DETERMINANT<T_ADAPTIVE,3>,typename ADAPTIVE_VECTOR_OPERATION_POLICY<NIL,T_ADAPTIVE,3>::DETERMINANT>
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
public:
    typedef typename ADAPTIVE_WRAPPER_BASE<ADAPTIVE_DETERMINANT,typename ADAPTIVE_VECTOR_OPERATION_POLICY<NIL,T_ADAPTIVE,3>::DETERMINANT>::WRAPPED_TYPE WRAPPED_TYPE;
private:
    const T_ADAPTIVE a00,a01,a02,a10,a11,a12,a20,a21,a22;
    ADAPTIVE_DETERMINANT();
    void operator=(const ADAPTIVE_DETERMINANT&);

public:
    ADAPTIVE_DETERMINANT(const T_ADAPTIVE& a11_input,const T_ADAPTIVE& a12_input,const T_ADAPTIVE& a13_input,const T_ADAPTIVE& a21_input,const T_ADAPTIVE& a22_input,
        const T_ADAPTIVE& a23_input,const T_ADAPTIVE& a31_input,const T_ADAPTIVE& a32_input,const T_ADAPTIVE& a33_input)
        :a00(a11_input),a01(a12_input),a02(a13_input),a10(a21_input),a11(a22_input),a12(a23_input),a20(a31_input),a21(a32_input),a22(a33_input)
    {}

    WRAPPED_TYPE Wrapped_Implementation() const
    {return (a00*(a11*a22-a12*a21)+a01*(a12*a20-a10*a22)+a02*(a10*a21-a11*a20));}
//#####################################################################
};

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Determinant(const T& t00)
{return ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>(t00);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Determinant(const VECTOR<T,1>& t1)
{return Adaptive_Determinant<T_EXACT>(t1[0]);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,1>
Adaptive_Determinant(const VECTOR<VECTOR<T,1>,1>& t)
{return Adaptive_Determinant<T_EXACT>(t[0]);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Determinant(const T& t00,const T& t01,const T& t10,const T& t11)
{return ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>(t00,t01,t10,t11);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Determinant(const VECTOR<T,2>& t1,const VECTOR<T,2>& t2)
{return Adaptive_Determinant<T_EXACT>(t1[0],t2[0],t1[1],t2[1]);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,2>
Adaptive_Determinant(const VECTOR<VECTOR<T,2>,2>& t)
{return Adaptive_Determinant<T_EXACT>(t[0],t[1]);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Determinant(const T& t00,const T& t01,const T& t02,const T& t10,const T& t11,const T& t12,const T& t20,const T& t21,const T& t22)
{return ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>(t00,t01,t02,t10,t11,t12,t20,t21,t22);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Determinant(const VECTOR<T,3>& t1,const VECTOR<T,3>& t2,const VECTOR<T,3>& t3)
{return Adaptive_Determinant<T_EXACT>(t1[0],t2[0],t3[0],t1[1],t2[1],t3[1],t1[2],t2[2],t3[2]);}

template<class T_EXACT,class T>
ADAPTIVE_DETERMINANT<typename ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE,3>
Adaptive_Determinant(const VECTOR<VECTOR<T,3>,3>& t)
{return Adaptive_Determinant<T_EXACT>(t[0],t[1],t[2]);}

}
using ADAPTIVE_DETAIL::ADAPTIVE_DETERMINANT;
using ADAPTIVE_DETAIL::Adaptive_Determinant;
}
#endif
