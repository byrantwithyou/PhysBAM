//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Igor Neverov, Craig Schroeder, Eftychios Sifakis, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUATERNION
//#####################################################################
#ifndef __QUATERNION__
#define __QUATERNION__

#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T> struct IS_SCALAR_BLOCK<QUATERNION<T> >:public IS_SCALAR_BLOCK<T>{};
template<class T> struct IS_SCALAR_VECTOR_SPACE<QUATERNION<T> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<QUATERNION<T>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class QUATERNION
{
    typedef VECTOR<T,3> TV;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T SCALAR;

    T s;
    TV v;

    QUATERNION()
        :s(0) // note that v is also zeroed
    {}

    template<class T2> explicit QUATERNION(const QUATERNION<T2>& q)
        :s((T)q.s),v(q.v)
    {}

    QUATERNION(const T s,const T x,const T y,const T z)
        :s(s),v(x,y,z)
    {}

    QUATERNION(const T s,const TV& v)
        :s(s),v(v)
    {}

    explicit QUATERNION(const VECTOR<T,4>& q)
        :s(q[0]),v(q[1],q[2],q[3])
    {}

    VECTOR<T,4> Vector() const
    {return VECTOR<T,4>(s,v.x,v.y,v.z);}

    static QUATERNION One()
    {return QUATERNION(1,0,0,0);}

    bool operator==(const QUATERNION& q) const
    {return s==q.s && v==q.v;}

    bool operator!=(const QUATERNION& q) const
    {return s!=q.s || v!=q.v;}

    QUATERNION operator-() const
    {return QUATERNION(-s,-v);}

    QUATERNION& operator+=(const QUATERNION& q)
    {s+=q.s;v+=q.v;return *this;}

    QUATERNION& operator-=(const QUATERNION& q)
    {s-=q.s;v-=q.v;return *this;}

    QUATERNION& operator*=(const QUATERNION& q)
    {return *this=*this*q;}

    QUATERNION& operator*=(const T a)
    {s*=a;v*=a;return *this;}

    QUATERNION& operator/=(const T a)
    {assert(a!=0);T r=1/a;s*=r;v*=r;return *this;}

    QUATERNION operator+(const QUATERNION& q) const
    {return QUATERNION(s+q.s,v+q.v);}

    QUATERNION operator-(const QUATERNION& q) const
    {return QUATERNION(s-q.s,v-q.v);}

    QUATERNION operator*(const QUATERNION& q) const // 16 mult and 13 add/sub
    {return QUATERNION(s*q.s-TV::Dot_Product(v,q.v),s*q.v+q.s*v+TV::Cross_Product(v,q.v));}

    QUATERNION operator*(const T a) const
    {return QUATERNION(s*a,v*a);}

    QUATERNION operator/(const T a) const
    {assert(a!=0);T r=1/a;return QUATERNION(s*r,v*r);}

    T Magnitude() const
    {return sqrt(Magnitude_Squared());}

    T Magnitude_Squared() const
    {return sqr(s)+v.Magnitude_Squared();}

    T Max_Abs() const
    {return maxabs(s,v.Max_Abs());}

    T Normalize()
    {T magnitude=Magnitude();if(magnitude) *this/=magnitude;else *this=One();return magnitude;}

    QUATERNION Normalized() const
    {QUATERNION q(*this);q.Normalize();return q;}

    bool Is_Normalized(const T tolerance=(T)1e-3) const
    {return abs(Magnitude_Squared()-(T)1)<=tolerance;}

    QUATERNION Conjugate() const
    {return QUATERNION(s,-v);}

    QUATERNION Inverse() const
    {return Conjugate()/Magnitude_Squared();}

    static T Dot_Product(const QUATERNION& q1,const QUATERNION& q2)
    {return q1.s*q2.s+TV::Dot_Product(q1.v,q2.v);}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,s,v);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,s,v);}

//#####################################################################
};

template<class T>
inline QUATERNION<T> operator*(const T a,const QUATERNION<T>& q)
{return QUATERNION<T>(q.s*a,q.v*a);}

template<class T>
inline std::ostream& operator<<(std::ostream& output,const QUATERNION<T>& q)
{output<<"("<<q.s<<" "<<q.v<<")";return output;}

template<class T>
inline std::istream& operator>>(std::istream& input,QUATERNION<T>& q)
{Ignore(input,'(');input>>q.s>>q.v;Ignore(input,')');return input;}

}
#endif
