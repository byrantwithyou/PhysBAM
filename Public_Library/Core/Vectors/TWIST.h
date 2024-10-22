//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TWIST
//#####################################################################
#ifndef __TWIST__
#define __TWIST__

#include <Core/Vectors/VECTOR.h>
#ifdef DIFFERENCE // Windows workaround.
#undef DIFFERENCE
#endif
namespace PhysBAM{

template<class TV> struct IS_SCALAR_BLOCK<TWIST<TV> >:public IS_SCALAR_BLOCK<TV>{using IS_SCALAR_BLOCK<TV>::value;};
template<class TV> struct IS_SCALAR_VECTOR_SPACE<TWIST<TV> >:public IS_SCALAR_VECTOR_SPACE<TV>{};
template<class TV,class RW> struct IS_BINARY_IO_SAFE<TWIST<TV>,RW>:public IS_BINARY_IO_SAFE<TV,RW>{};

template<class TV>
class TWIST
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
public:
    typedef T SCALAR;
    typedef int HAS_UNTYPED_READ_WRITE;

    enum WORKAROUND {dimension=TV::m+T_SPIN::m,m=dimension};

    TV linear;
    T_SPIN angular;
    
    TWIST()
        :angular()
    {}

    TWIST(const TV& linear_input,const T_SPIN& angular_input)
        :linear(linear_input),angular(angular_input)
    {}

    template<class T2> explicit TWIST(const TWIST<VECTOR<T2,TV::m> >& twist_input)
        :linear((TV)twist_input.linear),angular((T_SPIN)twist_input.angular)
    {}

    bool operator==(const TWIST& v) const
    {return linear==v.linear && angular==v.angular;}

    bool operator!=(const TWIST& v) const
    {return !(*this==v);}

    TWIST& operator+=(const TWIST& v)
    {linear+=v.linear;angular+=v.angular;return *this;}

    TWIST& operator-=(const TWIST& v)
    {linear-=v.linear;angular-=v.angular;return *this;}

    TWIST& operator*=(const T a)
    {linear*=a;angular*=a;return *this;}

    TWIST operator-() const
    {return TWIST(-linear,-angular);}

    TWIST operator+(const TWIST& v) const
    {return TWIST(linear+v.linear,angular+v.angular);}

    TWIST operator-(const TWIST& v) const
    {return TWIST(linear-v.linear,angular-v.angular);}

    TWIST operator*(const T a) const
    {return TWIST<TV>(linear*a,angular*a);}

    TWIST operator/(const T a) const
    {return *this*(1/a);}

    VECTOR<T,dimension> Get_Vector() const
    {return VECTOR<T,dimension>(linear,angular);}

    void Set_Vector(const VECTOR<T,dimension>& vector)
    {vector.Extract(linear,angular);}

    T Dot(const TWIST& v) const
    {return TV::Dot_Product(linear,v.linear)+T_SPIN::Dot_Product(angular,v.angular);};

    static T Dot_Product(const TWIST& v1,const TWIST& v2)
    {return v1.Dot(v2);}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,linear,angular);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,linear,angular);}
//#####################################################################
};
// global functions
template<class TV> inline TWIST<TV> operator*(const typename TV::SCALAR a,const TWIST<TV>& v)
{return TWIST<TV>(a*v.linear,a*v.angular);}

//#####################################################################
template<class TV> struct CAN_ASSIGN<TWIST<TV>,TWIST<TV> > {static const bool value=true;};
//#####################################################################

template<class TV> inline std::istream& operator>>(std::istream& input,TWIST<TV>& v)
{Ignore(input,'(');input>>v.linear>>v.angular;Ignore(input,')');return input;}

template<class TV> inline std::ostream& operator<<(std::ostream& output,const TWIST<TV>& v)
{output<<"("<<v.linear<<"  "<<v.angular<<")";return output;}
//#####################################################################
// Function Random_Fill_Uniform
//#####################################################################
template<class T,class TV> void
Random_Fill_Uniform(RANDOM_NUMBERS<T>& rand,TWIST<TV>& m,const T a,const T b)
{
    rand.Fill_Uniform(m.linear,a,b);
    rand.Fill_Uniform(m.angular,a,b);
}
}
#endif
