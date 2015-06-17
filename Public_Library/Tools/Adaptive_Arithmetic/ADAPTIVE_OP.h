//#####################################################################
// Copyright 2008, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_OP__
#define __ADAPTIVE_OP__
#include <Tools/Adaptive_Arithmetic/ADAPTIVE_ATOM.h>
#include <Tools/Adaptive_Arithmetic/ADAPTIVE_BASE.h>
#include <Tools/Adaptive_Arithmetic/ADAPTIVE_DIFFERENCE.h>
#include <Tools/Adaptive_Arithmetic/ADAPTIVE_NEGATION.h>
#include <Tools/Adaptive_Arithmetic/ADAPTIVE_PRODUCT.h>
#include <Tools/Adaptive_Arithmetic/ADAPTIVE_QUOTIENT.h>
#include <Tools/Adaptive_Arithmetic/ADAPTIVE_SUM.h>
#include <Tools/Adaptive_Arithmetic/ADAPTIVE_VECTOR_OPERATION_POLICY.h>
#include <Tools/Adaptive_Arithmetic/EXACT_ARITHMETIC_POLICY.h>
#include <Tools/Math_Tools/ZERO.h>
#include <Tools/Utilities/STATIC_ASSERT.h>
#include <Tools/Vectors/VECTOR_3D.h>
#include <iostream>

namespace PhysBAM{
//#####################################################################
// Adaptive_Vector_Sum
//#####################################################################
template<class T_EXACT,class T> typename ADAPTIVE_VECTOR_OPERATION_POLICY<T_EXACT,T,2>::SUM
Adaptive_Vector_Sum(const VECTOR<T,2>& v)
{
    typedef typename ADAPTIVE_DETAIL::ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE T_ADAPTIVE;
    return T_ADAPTIVE(v[0])+T_ADAPTIVE(v[1]);
}
template<class T_EXACT,class T> typename ADAPTIVE_VECTOR_OPERATION_POLICY<T_EXACT,T,3>::SUM
Adaptive_Vector_Sum(const VECTOR<T,3>& v)
{
    typedef typename ADAPTIVE_DETAIL::ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE T_ADAPTIVE;
    return T_ADAPTIVE(v[0])+T_ADAPTIVE(v[1])+T_ADAPTIVE(v[2]);
}
template<class T_EXACT,class T> typename ADAPTIVE_VECTOR_OPERATION_POLICY<T_EXACT,T,4>::SUM
Adaptive_Vector_Sum(const VECTOR<T,4>& v)
{
    typedef typename ADAPTIVE_DETAIL::ADAPTIVE_POLICY<T_EXACT,T>::ADAPTIVE T_ADAPTIVE;
    return T_ADAPTIVE(v[0])+T_ADAPTIVE(v[1])+T_ADAPTIVE(v[2])+T_ADAPTIVE(v[3]);
}
//#####################################################################
// operator-
//#####################################################################
template<class T_ADAPTIVE0>
ADAPTIVE_NEGATION<T_ADAPTIVE0>
operator-(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return ADAPTIVE_NEGATION<T_ADAPTIVE0>(a1.Derived());
}
//#####################################################################
// operator+
//#####################################################################
template<class T_ADAPTIVE0,class T_ADAPTIVE1>
ADAPTIVE_SUM<T_ADAPTIVE0,T_ADAPTIVE1>
operator+(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const ADAPTIVE_BASE<T_ADAPTIVE1,typename T_ADAPTIVE1::EXACT_TYPE>& a2)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE1>::value));
    return ADAPTIVE_SUM<T_ADAPTIVE0,T_ADAPTIVE1>(a1.Derived(),a2.Derived());
}

template<class T_ADAPTIVE0,class T>
typename enable_if<IS_ATOMIZABLE<T,typename T_ADAPTIVE0::EXACT_TYPE>::value,ADAPTIVE_SUM<T_ADAPTIVE0,ADAPTIVE_ATOM<T,typename T_ADAPTIVE0::EXACT_TYPE> > >::type
operator+(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const T& t)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    typedef typename T_ADAPTIVE0::EXACT_TYPE EXACT_TYPE;
    STATIC_ASSERT((IS_ATOMIZABLE<T,EXACT_TYPE>::value));
    typedef ADAPTIVE_ATOM<T,EXACT_TYPE> ATOM;
    return ADAPTIVE_SUM<T_ADAPTIVE0,ATOM>(a1.Derived(),ATOM(t));
}

template<class T,class T_ADAPTIVE0>
typename enable_if<IS_ATOMIZABLE<T,typename T_ADAPTIVE0::EXACT_TYPE>::value,ADAPTIVE_SUM<ADAPTIVE_ATOM<T,typename T_ADAPTIVE0::EXACT_TYPE>,T_ADAPTIVE0> >::type
operator+(const T& t,const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    typedef typename T_ADAPTIVE0::EXACT_TYPE EXACT_TYPE;
    STATIC_ASSERT((IS_ATOMIZABLE<T,EXACT_TYPE>::value));
    typedef ADAPTIVE_ATOM<T,EXACT_TYPE> ATOM;
    return ADAPTIVE_SUM<ATOM,T_ADAPTIVE0>(ATOM(t),a1.Derived());
}

//#####################################################################
// operator-
//#####################################################################
template<class T_ADAPTIVE0,class T_ADAPTIVE1>
ADAPTIVE_DIFFERENCE<T_ADAPTIVE0,T_ADAPTIVE1>
operator-(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const ADAPTIVE_BASE<T_ADAPTIVE1,typename T_ADAPTIVE1::EXACT_TYPE>& a2)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE1>::value));
    return ADAPTIVE_DIFFERENCE<T_ADAPTIVE0,T_ADAPTIVE1>(a1.Derived(),a2.Derived());
}

template<class T_ADAPTIVE0,class T>
typename enable_if<IS_ATOMIZABLE<T,typename T_ADAPTIVE0::EXACT_TYPE>::value,
ADAPTIVE_DIFFERENCE<T_ADAPTIVE0,ADAPTIVE_ATOM<T,typename T_ADAPTIVE0::EXACT_TYPE> > >::TYPE
operator-(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const T& t)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    typedef typename T_ADAPTIVE0::EXACT_TYPE EXACT_TYPE;
    STATIC_ASSERT((IS_ATOMIZABLE<T,EXACT_TYPE>::value));
    typedef ADAPTIVE_ATOM<T,EXACT_TYPE> ATOM;
    return ADAPTIVE_DIFFERENCE<T_ADAPTIVE0,ATOM>(a1.Derived(),ATOM(t));
}

template<class T,class T_ADAPTIVE0>
typename enable_if<IS_ATOMIZABLE<T,typename T_ADAPTIVE0::EXACT_TYPE>::value,
ADAPTIVE_DIFFERENCE<ADAPTIVE_ATOM<T,typename T_ADAPTIVE0::EXACT_TYPE>,T_ADAPTIVE0> >::TYPE
operator-(const T& t,const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    typedef typename T_ADAPTIVE0::EXACT_TYPE EXACT_TYPE;
    STATIC_ASSERT((IS_ATOMIZABLE<T,EXACT_TYPE>::value));
    typedef ADAPTIVE_ATOM<T,EXACT_TYPE> ATOM;
    return ADAPTIVE_DIFFERENCE<ATOM,T_ADAPTIVE0>(ATOM(t),a1.Derived());
}

//#####################################################################
// operator*
//#####################################################################
template<class T_ADAPTIVE0,class T_ADAPTIVE1>
ADAPTIVE_PRODUCT<T_ADAPTIVE0,T_ADAPTIVE1>
operator*(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const ADAPTIVE_BASE<T_ADAPTIVE1,typename T_ADAPTIVE1::EXACT_TYPE>& a2)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE1>::value));
    return ADAPTIVE_PRODUCT<T_ADAPTIVE0,T_ADAPTIVE1>(a1.Derived(),a2.Derived());
}

template<class T_ADAPTIVE0,class T>
typename enable_if<IS_ATOMIZABLE<T,typename T_ADAPTIVE0::EXACT_TYPE>::value,
ADAPTIVE_PRODUCT<T_ADAPTIVE0,ADAPTIVE_ATOM<T,typename T_ADAPTIVE0::EXACT_TYPE> > >::TYPE
operator*(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const T& t)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    typedef typename T_ADAPTIVE0::EXACT_TYPE EXACT_TYPE;
    STATIC_ASSERT((IS_ATOMIZABLE<T,EXACT_TYPE>::value));
    typedef ADAPTIVE_ATOM<T,EXACT_TYPE> ATOM;
    return ADAPTIVE_PRODUCT<T_ADAPTIVE0,ATOM>(a1.Derived(),ATOM(t));
}

template<class T,class T_ADAPTIVE0>
typename enable_if<IS_ATOMIZABLE<T,typename T_ADAPTIVE0::EXACT_TYPE>::value,
ADAPTIVE_PRODUCT<ADAPTIVE_ATOM<T,typename T_ADAPTIVE0::EXACT_TYPE>,T_ADAPTIVE0> >::TYPE
operator*(const T& t,const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    typedef typename T_ADAPTIVE0::EXACT_TYPE EXACT_TYPE;
    STATIC_ASSERT((IS_ATOMIZABLE<T,EXACT_TYPE>::value));
    typedef ADAPTIVE_ATOM<T,EXACT_TYPE> ATOM;
    return ADAPTIVE_PRODUCT<ATOM,T_ADAPTIVE0>(ATOM(t),a1.Derived());
}


//#####################################################################
// operator/
//#####################################################################
template<class T_ADAPTIVE0,class T_ADAPTIVE1>
ADAPTIVE_QUOTIENT<T_ADAPTIVE0,T_ADAPTIVE1>
operator/(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const ADAPTIVE_BASE<T_ADAPTIVE1,typename T_ADAPTIVE1::EXACT_TYPE>& a2)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE1>::value));
    return ADAPTIVE_QUOTIENT<T_ADAPTIVE0,T_ADAPTIVE1>(a1.Derived(),a2.Derived());
}

template<class T_ADAPTIVE0,class T>
typename enable_if<IS_ATOMIZABLE<T,typename T_ADAPTIVE0::EXACT_TYPE>::value,
ADAPTIVE_QUOTIENT<T_ADAPTIVE0,ADAPTIVE_ATOM<T,typename T_ADAPTIVE0::EXACT_TYPE> > >::TYPE
operator/(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const T& t)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    typedef typename T_ADAPTIVE0::EXACT_TYPE EXACT_TYPE;
    STATIC_ASSERT((IS_ATOMIZABLE<T,EXACT_TYPE>::value));
    typedef ADAPTIVE_ATOM<T,EXACT_TYPE> ATOM;
    return ADAPTIVE_QUOTIENT<T_ADAPTIVE0,ATOM>(a1.Derived(),ATOM(t));
}

template<class T,class T_ADAPTIVE0>
typename enable_if<IS_ATOMIZABLE<T,typename T_ADAPTIVE0::EXACT_TYPE>::value,
ADAPTIVE_QUOTIENT<ADAPTIVE_ATOM<T,typename T_ADAPTIVE0::EXACT_TYPE>,T_ADAPTIVE0> >::TYPE
operator/(const T& t,const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    typedef typename T_ADAPTIVE0::EXACT_TYPE EXACT_TYPE;
    STATIC_ASSERT((IS_ATOMIZABLE<T,EXACT_TYPE>::value));
    typedef ADAPTIVE_ATOM<T,EXACT_TYPE> ATOM;
    return ADAPTIVE_QUOTIENT<ATOM,T_ADAPTIVE0>(ATOM(t),a1.Derived());
}

//#####################################################################
// Compare
//#####################################################################
template<class T_ADAPTIVE0>
int Compare(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,ZERO zero)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return a1.Sign();
}

template<class T_ADAPTIVE0>
int Compare(ZERO zero,const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return -Compare(a1,zero);
}

template<class T_ADAPTIVE0,class T>
typename enable_if<IS_ATOMIZABLE<T,typename T_ADAPTIVE0::EXACT_TYPE>::value,int>::type
Compare(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const T& t)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return Compare(a1,ADAPTIVE_ATOM<T,typename T_ADAPTIVE0::EXACT_TYPE>(t));
}

template<class T,class T_ADAPTIVE0>
typename enable_if<IS_ATOMIZABLE<T,typename T_ADAPTIVE0::EXACT_TYPE>::value,int>::type
Compare(const T& t,const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return -Compare(t,a1);
}

template<class T_ADAPTIVE0,class T_ADAPTIVE1>
int Compare(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const ADAPTIVE_BASE<T_ADAPTIVE1,typename T_ADAPTIVE1::EXACT_TYPE>& a2)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE1>::value));
    return Compare(a1-a2,ZERO());
}

//#####################################################################
// operator<
//#####################################################################
template<class T_ADAPTIVE0,class T> 
bool operator<(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const T& t)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return (Compare(a1,t)<0);
}

template<class T,class T_ADAPTIVE0>
typename enable_if<IS_NOT_ADAPTIVE<T>::value,bool>::type operator<(const T& t,const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return (Compare(t,a1)<0);
}

//#####################################################################
// operator<=
//#####################################################################
template<class T_ADAPTIVE0,class T>
bool operator<=(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const T& t)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return !operator<(t,a1);
}

template<class T,class T_ADAPTIVE0>
typename enable_if<IS_NOT_ADAPTIVE<T>::value,bool>::type 
operator<=(const T& t,const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return !operator<(a1,t);
}

//#####################################################################
// operator>=
//#####################################################################
template<class T_ADAPTIVE0,class T>
bool operator>=(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const T& t)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return !operator<(a1,t);
}

template<class T,class T_ADAPTIVE0>
typename enable_if<IS_NOT_ADAPTIVE<T>::value,
bool>::TYPE
operator>=(const T& t,const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return !operator<(t,a1);
}

//#####################################################################
// operator>
//#####################################################################
template<class T_ADAPTIVE0,class T>
bool 
operator>(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const T& t)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return operator<(t,a1);
}

template<class T,class T_ADAPTIVE0>
typename enable_if<IS_NOT_ADAPTIVE<T>::value,bool>::type
operator>(const T& t,const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return operator<(a1,t);
}

//#####################################################################
// operator==
//#####################################################################
template<class T_ADAPTIVE0,class T>
bool operator==(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const T& t)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return (Compare(a1,t)==0);
}

template<class T,class T_ADAPTIVE0>
typename enable_if<IS_NOT_ADAPTIVE<T>::value,bool>::type
operator==(const T& t,const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return operator==(a1,t);
}

//#####################################################################
// operator!=
//#####################################################################
template<class T_ADAPTIVE0,class T>
bool operator!=(const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1,const T& t)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return !operator==(a1,t);
}

template<class T,class T_ADAPTIVE0>
typename enable_if<IS_NOT_ADAPTIVE<T>::value,bool>::type
operator!=(const T& t,const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    return !operator==(t,a1);
}

//#####################################################################
// operator<<
//#####################################################################
template<class T_ADAPTIVE0>
std::ostream&
operator<<(std::ostream& out,const ADAPTIVE_BASE<T_ADAPTIVE0,typename T_ADAPTIVE0::EXACT_TYPE>& a1)
{
    typedef typename T_ADAPTIVE0::FP_TYPE FP_TYPE;
    FP_TYPE estimate,error;
    a1.Estimate_And_Error().Get(estimate,error);
    return out <<"(" <<estimate <<"+/- " <<error <<")";
}
//#####################################################################
}
#endif
