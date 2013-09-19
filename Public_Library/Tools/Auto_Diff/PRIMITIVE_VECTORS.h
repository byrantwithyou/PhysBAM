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

template<class T_OBJ,class ENABLER=void> struct GET_FLAGS;
template<class T_OBJ> struct GET_FLAGS<T_OBJ,typename IF<1,void,typename T_OBJ::FLAGS>::TYPE> {typedef typename T_OBJ::FLAGS TYPE;};

template<int m> struct V_FLAGS;

template<class TV>
struct ZERO_VEC
{
    typedef V_FLAGS<0> FLAGS;
    typedef typename TV::SCALAR T;
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

template<int m> struct V_FLAGS
{
    enum WA{mask=m};
    V_FLAGS<m> operator+() const;
    V_FLAGS<m> operator-() const;
    template<int p> V_FLAGS<m|p> operator+(V_FLAGS<p> x) const;
    template<int p> V_FLAGS<m|p> operator-(V_FLAGS<p> x) const;
    V_FLAGS<m> operator*(double a) const;
    V_FLAGS<m> operator/(double a) const;
};
template<class T,int d> struct GET_FLAGS<VECTOR<T,d> > {typedef V_FLAGS<1> TYPE;};
template<int m,int m2> V_FLAGS<m|m2> Choose(V_FLAGS<m>,V_FLAGS<m2>);

template<class TV> TV operator+ (const ZERO_VEC<TV>& z,const TV& v) {return v;}
template<class TV> TV operator+ (const TV& v,const ZERO_VEC<TV>& z) {return v;}
template<class TV> TV operator- (const ZERO_VEC<TV>& z,const TV& v) {return -v;}
template<class TV> TV operator- (const TV& v,const ZERO_VEC<TV>& z) {return v;}

template<class TV,class G> struct CH_VEC;
template<class TV> struct CH_VEC<TV,V_FLAGS<0> > {typedef ZERO_VEC<TV> TYPE;};
template<class TV> struct CH_VEC<TV,V_FLAGS<1> > {typedef TV TYPE;};
#define MK_VEC(TV,G) typename CH_VEC<TV,decltype(G)>::TYPE

template<class T> void Fill_From_Helper(T& a,const T& b){a=b;}
template<class T,int d> void Fill_From_Helper(VECTOR<T,d>& m,const ZERO_VEC<VECTOR<T,d> >& z){m=VECTOR<T,d>();}
template<class T,int d> inline VECTOR<T,d> Cast_Helper(const VECTOR<T,d>& z){return z;}
template<class TV> inline TV Cast_Helper(ZERO_VEC<TV> z){return TV();}
}
}
#endif
