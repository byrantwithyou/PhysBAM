//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PRIMITIVE_VECTORS
//##################################################################### 
#ifndef __PRIMITIVE_VECTORS__
#define __PRIMITIVE_VECTORS__

#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{
namespace HETERO_DIFF{

template<class TV>
struct ZERO_VEC
{
    typedef typename TV::SCALAR SCALAR;typedef SCALAR T;
    ZERO_VEC operator-() const
    {return *this;}

    ZERO_VEC operator+(const ZERO_VEC&) const
    {return *this;}

    ZERO_VEC operator-(const ZERO_VEC&) const
    {return *this;}

    ZERO_VEC operator*(const T&) const
    {return *this;}

    ZERO_VEC operator/(const T&) const
    {return *this;}

    T Dot(const ZERO_VEC&) const
    {return T();}
};

template<class T> struct IS_VECTOR{static const int value=0;};
template<class TV> struct IS_VECTOR<ZERO_VEC<TV> > {static const int value=1;};
template<class T,int d> struct IS_VECTOR<VECTOR<T,d> > {static const int value=1;};

template<class TV> TV operator+ (const ZERO_VEC<TV>& z,const TV& v) {return v;}
template<class TV> TV operator+ (const TV& v,const ZERO_VEC<TV>& z) {return v;}
template<class TV> TV operator- (const ZERO_VEC<TV>& z,const TV& v) {return -v;}
template<class TV> TV operator- (const TV& v,const ZERO_VEC<TV>& z) {return v;}

template<class T> typename ENABLE_IF<IS_VECTOR<T>::value>::TYPE Fill_From(T& a,const T& b){a=b;}
template<class T,int d> void Fill_From(VECTOR<T,d>& m,const ZERO_VEC<VECTOR<T,d> >& z){m=VECTOR<T,d>();}
template<class T,int d> inline VECTOR<T,d> Cast_Helper(const VECTOR<T,d>& z){return z;}
template<class TV> inline TV Cast_Helper(ZERO_VEC<TV> z){return TV();}

template<class T_VEC> typename ENABLE_IF<IS_VECTOR<T_VEC>::value,T_VEC>::TYPE Choose(const T_VEC& a,const T_VEC& b);
template<class TV,class T_VEC> T_VEC Choose(const T_VEC& a,const ZERO_VEC<TV>& b);
template<class TV,class T_VEC> T_VEC Choose(const ZERO_VEC<TV>& a,const T_VEC& b);
template<class TV> ZERO_VEC<TV> Choose(const ZERO_VEC<TV>& a,const ZERO_VEC<TV>& b);
}
}
#endif
