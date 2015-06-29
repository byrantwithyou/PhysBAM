//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ZERO_VECTOR
//##################################################################### 
#ifndef __ZERO_VECTOR__
#define __ZERO_VECTOR__

#include <Tools/Math_Tools/FIXED_NUMBER.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

template<class T,int d>
class ZERO_VECTOR
{
public:
    enum {m=d};
    typedef T SCALAR;
    ZERO_VECTOR operator+() const
    {return *this;}

    ZERO_VECTOR operator-() const
    {return *this;}

    ZERO_VECTOR operator+(const ZERO_VECTOR&) const
    {return *this;}

    ZERO_VECTOR operator-(const ZERO_VECTOR&) const
    {return *this;}

    ZERO_VECTOR operator*(const T&) const
    {return *this;}

    ZERO_VECTOR operator/(const T&) const
    {return *this;}

    FIXED_NUMBER<T,0> Dot(const ZERO_VECTOR&) const
    {return FIXED_NUMBER<T,0>();}

    FIXED_NUMBER<T,0> Dot(const VECTOR<T,d>&) const
    {return FIXED_NUMBER<T,0>();}
};

template<class T,int d> VECTOR<T,d> operator+=(const VECTOR<T,d>& v,const ZERO_VECTOR<T,d>& z){return v;}

template<class T,int d> VECTOR<T,d> operator+ (const ZERO_VECTOR<T,d>& z,const VECTOR<T,d>& v) {return v;}
template<class T,int d> VECTOR<T,d> operator+ (const VECTOR<T,d>& v,const ZERO_VECTOR<T,d>& z) {return v;}
template<class T,int d> VECTOR<T,d> operator- (const ZERO_VECTOR<T,d>& z,const VECTOR<T,d>& v) {return -v;}
template<class T,int d> VECTOR<T,d> operator- (const VECTOR<T,d>& v,const ZERO_VECTOR<T,d>& z) {return v;}

template<class T,int d> void Fill_From(VECTOR<T,d>& m,const ZERO_VECTOR<T,d>& z){m=VECTOR<T,d>();}
template<class T,int d> inline VECTOR<T,d> Cast_Helper(const VECTOR<T,d>& z){return z;}
template<class T,int d> inline VECTOR<T,d> Cast_Helper(ZERO_VECTOR<T,d> z){return VECTOR<T,d>();}

template<class T_VEC> typename enable_if<IS_VECTOR<T_VEC>::value,T_VEC>::type Choose(const T_VEC& a,const T_VEC& b);
template<class T,int d,class T_VEC> T_VEC Choose(const T_VEC& a,const ZERO_VECTOR<T,d>& b);
template<class T,int d,class T_VEC> T_VEC Choose(const ZERO_VECTOR<T,d>& a,const T_VEC& b);
template<class T,int d> ZERO_VECTOR<T,d> Choose(const ZERO_VECTOR<T,d>& a,const ZERO_VECTOR<T,d>& b);

template<class T,int d> inline ZERO_VECTOR<T,d> operator*(const FIXED_NUMBER<T,0>,const VECTOR<T,d>&) {return ZERO_VECTOR<T,d>();}
template<class T,int d> inline ZERO_VECTOR<T,d> operator/(const FIXED_NUMBER<T,0>,const VECTOR<T,d>&) {return ZERO_VECTOR<T,d>();}
template<class T,int d> inline ZERO_VECTOR<T,d> operator*(const VECTOR<T,d>&,const FIXED_NUMBER<T,0>) {return ZERO_VECTOR<T,d>();}

template<class T,int d> inline typename enable_if<is_scalar<T>::value,ZERO_VECTOR<T,d> >::type operator*(const T&,const ZERO_VECTOR<T,d>&) {return ZERO_VECTOR<T,d>();}

template<class T> ZERO_VECTOR<T,3> Cross_Product(const ZERO_VECTOR<T,3>&,const ZERO_VECTOR<T,3>&){return ZERO_VECTOR<T,3>();}
template<class T> ZERO_VECTOR<T,3> Cross_Product(const ZERO_VECTOR<T,3>&,const VECTOR<T,3>&){return ZERO_VECTOR<T,3>();}
template<class T> ZERO_VECTOR<T,3> Cross_Product(const VECTOR<T,3>&,const ZERO_VECTOR<T,3>&){return ZERO_VECTOR<T,3>();}
template<class T> VECTOR<T,3> Cross_Product(const VECTOR<T,3>& a,const VECTOR<T,3>& b){return a.Cross(b);}
}
#endif
