//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUANTITY
//#####################################################################
#ifndef __QUANTITY__
#define __QUANTITY__

#include <PhysBAM_Tools/Math_Tools/cbrt.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <cfloat>
#include <cmath>
#include <iosfwd>
#include "CASTER.h"
#include "PARTIAL_UNIT.h"
#include "QUANTITY_FORWARD.h"
#include <boost/math/special_functions/asinh.hpp>
namespace PhysBAM{
namespace UNITS{

template<class T>
class QUANTITY
{
    struct UNUSABLE{};
public:
    T value;
    PARTIAL_UNIT unit; 

    QUANTITY()
        :value(),unit(PARTIAL_UNIT::Zero())
    {}

    QUANTITY(const int a)
        :value(a),unit(PARTIAL_UNIT::One())
    {}

    QUANTITY(const unsigned a)
        :value(a),unit(PARTIAL_UNIT::One())
    {}

    QUANTITY(const long a)
        :value(a),unit(PARTIAL_UNIT::One())
    {}

    explicit QUANTITY(const long long a)
        :value(a),unit(PARTIAL_UNIT::One())
    {}

    QUANTITY(const float a)
        :value(a),unit(PARTIAL_UNIT::New())
    {}

    QUANTITY(const double a)
        :value(a),unit(PARTIAL_UNIT::New())
    {}

    template<class T2>
    QUANTITY(const QUANTITY<T2>& a)
        :value(a.value),unit(a.unit)
    {}

    QUANTITY(const T value,const PARTIAL_UNIT& unit)
        :value(value),unit(unit)
    {}

private:
    struct SAFE_BOOL_HELPER{void F(){}};
    typedef void (SAFE_BOOL_HELPER::*SAFE_BOOL)();
public:

    operator SAFE_BOOL() const // allow conversion to bool without allowing conversion to int
    {return value?&SAFE_BOOL_HELPER::F:0;}

    bool operator==(const QUANTITY& a) const
    {Unify(unit,a.unit);return value==a.value;}

    bool operator==(const int a) const
    {Unify_One(unit);return value==a;}

    bool operator==(const float a) const
    {Unify_One(unit);return value==a;}

    bool operator==(const double a) const
    {Unify_One(unit);return value==a;}

    bool operator!=(const QUANTITY& a) const
    {Unify(unit,a.unit);return value!=a.value;}

    bool operator!=(const int a) const
    {Unify_One(unit);return value!=a;}

    bool operator!=(const float a) const
    {Unify_One(unit);return value!=a;}

    bool operator!=(const double a) const
    {Unify_One(unit);return value!=a;}

    bool operator<(const QUANTITY& a) const
    {Unify(unit,a.unit);return value<a.value;}

    friend bool operator<(const int a,const QUANTITY& b)
    {Unify_One(b.unit);return a<b.value;}

    friend bool operator<(const float a,const QUANTITY& b)
    {Unify_One(b.unit);return a<b.value;}

    friend bool operator<(const double a,const QUANTITY& b)
    {Unify_One(b.unit);return a<b.value;}

    bool operator>(const QUANTITY& a) const
    {Unify(unit,a.unit);return value>a.value;}

    friend bool operator>(const int a,const QUANTITY& b)
    {Unify_One(b.unit);return a>b.value;}

    friend bool operator>(const float a,const QUANTITY& b)
    {Unify_One(b.unit);return a>b.value;}

    friend bool operator>(const double a,const QUANTITY& b)
    {Unify_One(b.unit);return a>b.value;}

    bool operator<=(const QUANTITY& a) const
    {Unify(unit,a.unit);return value<=a.value;}

    friend bool operator<=(const int a,const QUANTITY& b)
    {Unify_One(b.unit);return a<=b.value;}

    friend bool operator<=(const float a,const QUANTITY& b)
    {Unify_One(b.unit);return a<=b.value;}

    friend bool operator<=(const double a,const QUANTITY& b)
    {Unify_One(b.unit);return a<=b.value;}

    bool operator>=(const QUANTITY& a) const
    {Unify(unit,a.unit);return value>=a.value;}

    QUANTITY operator+() const
    {return QUANTITY(+value,unit);}

    QUANTITY operator-() const
    {return QUANTITY(-value,unit);}

    QUANTITY operator+(const QUANTITY& a) const
    {Unify(unit,a.unit);return QUANTITY(value+a.value,unit);}

    QUANTITY& operator+=(const QUANTITY& a)
    {Unify(unit,a.unit);value+=a.value;return *this;}

    QUANTITY operator-(const QUANTITY& a) const
    {Unify(unit,a.unit);return QUANTITY(value-a.value,unit);}

    QUANTITY& operator-=(const QUANTITY& a)
    {Unify(unit,a.unit);value-=a.value;return *this;}

    QUANTITY operator*(const QUANTITY& a) const
    {return QUANTITY(value*a.value,unit*a.unit);}

    QUANTITY operator*(const int a) const
    {return QUANTITY(value*a,unit);}

    QUANTITY operator*(const long a) const
    {return QUANTITY(value*a,unit);}

    QUANTITY operator*(const float a) const
    {return QUANTITY(value*a,unit);}

    QUANTITY operator*(const double a) const
    {return QUANTITY(value*a,unit);}

    QUANTITY& operator*=(const QUANTITY& a)
    {unit*=a.unit;value*=a.value;return *this;}

    QUANTITY operator/(const QUANTITY& a) const
    {return QUANTITY(value/a.value,unit/a.unit);}

    QUANTITY& operator/=(const QUANTITY& a)
    {unit/=a.unit;value/=a.value;return *this;}

    operator STREAM_TYPE() const
    {return STREAM_TYPE(T());}

    template<class RW>
    void Read(std::istream& input)
    {Read_Binary<RW>(input,value);unit=PARTIAL_UNIT::New();}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Binary<RW>(output,value);}

//#####################################################################
};

//#####################################################################
// Free operators
//#####################################################################
template<class T> inline QUANTITY<T>
operator+(const int a,const QUANTITY<T>& b)
{Unify_One(b.unit);return QUANTITY<T>(a+b.value,b.unit);}

template<class T> inline QUANTITY<T>
operator+(const T a,const QUANTITY<T>& b)
{Unify_One(b.unit);return QUANTITY<T>(a+b.value,b.unit);}

template<class T> inline QUANTITY<T>
operator-(const int a,const QUANTITY<T>& b)
{Unify_One(b.unit);return QUANTITY<T>(a-b.value,b.unit);}

template<class T> inline QUANTITY<T>
operator-(const float a,const QUANTITY<T>& b)
{Unify_One(b.unit);return QUANTITY<T>(a-b.value,b.unit);}

template<class T> inline QUANTITY<double>
operator-(const double a,const QUANTITY<T>& b)
{Unify_One(b.unit);return QUANTITY<T>(a-b.value,b.unit);}

template<class T> inline QUANTITY<T>
operator*(const int a,const QUANTITY<T>& b)
{return QUANTITY<T>(a*b.value,b.unit);}

template<class T> inline QUANTITY<T>
operator*(const float a,const QUANTITY<T>& b)
{return QUANTITY<T>(a*b.value,b.unit);}

template<class T> inline QUANTITY<double>
operator*(const double a,const QUANTITY<T>& b)
{return QUANTITY<T>(a*b.value,b.unit);}

template<class T> inline QUANTITY<T>
operator/(const int a,const QUANTITY<T>& b)
{return QUANTITY<T>(a/b.value,b.unit.Inverse());}

template<class T> inline QUANTITY<T>
operator/(const float a,const QUANTITY<T>& b)
{return QUANTITY<T>(a/b.value,b.unit.Inverse());}

template<class T> inline QUANTITY<double>
operator/(const double a,const QUANTITY<T>& b)
{return QUANTITY<T>(a/b.value,b.unit.Inverse());}

//#####################################################################
// Free functions
//#####################################################################
using ::std::abs;
using ::std::floor;
using ::std::ceil;
using ::std::sqrt;
using ::std::fmod;
using ::std::exp;
using ::std::log;
using ::std::pow;
using ::std::cos;
using ::std::sin;
using ::std::tan;
using ::std::acos;
using ::std::asin;
using ::std::atan;
using ::std::atan2;
using ::boost::math::asinh;

template<class T> inline QUANTITY<T>
abs(const QUANTITY<T>& a)
{return QUANTITY<T>(abs(a.value),a.unit);}

template<class T> inline QUANTITY<T>
floor(const QUANTITY<T>& a)
{Unify_One(a.unit);return QUANTITY<T>(floor(a.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
ceil(const QUANTITY<T>& a)
{Unify_One(a.unit);return QUANTITY<T>(ceil(a.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
sqrt(const QUANTITY<T>& a)
{return QUANTITY<T>(sqrt(a.value),sqrt(a.unit));}

inline QUANTITY<double>
cbrt(const QUANTITY<double> a)
{return QUANTITY<double>(PhysBAM::cbrt(a.value),cbrt(a.unit));}

template<class T> inline QUANTITY<T>
fmod(const QUANTITY<T>& a,const QUANTITY<T>& b)
{Unify(a.unit,b.unit);return QUANTITY<T>(fmod(a.value,b.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
exp(const QUANTITY<T>& a)
{Unify_One(a.unit);return QUANTITY<T>(exp(a.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
log(const QUANTITY<T>& a)
{Unify_One(a.unit);return QUANTITY<T>(log(a.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
cos(const QUANTITY<T>& a)
{Unify_One(a.unit);return QUANTITY<T>(cos(a.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
sin(const QUANTITY<T>& a)
{Unify_One(a.unit);return QUANTITY<T>(sin(a.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
tan(const QUANTITY<T>& a)
{Unify_One(a.unit);return QUANTITY<T>(tan(a.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
acos(const QUANTITY<T>& a)
{Unify_One(a.unit);return QUANTITY<T>(acos(a.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
asin(const QUANTITY<T>& a)
{Unify_One(a.unit);return QUANTITY<T>(asin(a.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
atan(const QUANTITY<T>& a)
{Unify_One(a.unit);return QUANTITY<T>(atan(a.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
asinh(const QUANTITY<T>& a)
{Unify_One(a.unit);return QUANTITY<T>(asinh(a.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
atan2(const QUANTITY<T>& x,const QUANTITY<T>& y)
{Unify(x.unit,y.unit);return QUANTITY<T>(atan2(x.value,y.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
pow(const QUANTITY<T>& x,const int y)
{return QUANTITY<T>(pow(x.value,y),pow(x.unit,y));}

template<class T> inline QUANTITY<double>
pow(const QUANTITY<T>& x,const double y)
{Unify_One(x.unit);return QUANTITY<T>(pow(x.value,y),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<double>
pow(const double x,const QUANTITY<T>& y)
{Unify_One(y.unit);return QUANTITY<T>(pow(x,(double)y.value),PARTIAL_UNIT::One());}

template<class T> inline QUANTITY<T>
pow(const QUANTITY<T>& x,const QUANTITY<T>& y)
{Unify_One(x.unit);Unify_One(y.unit);return QUANTITY<T>(pow(x.value,y.value),PARTIAL_UNIT::One());}

}

//#####################################################################
// Input and output
//#####################################################################
template<> inline void
Read_Binary<QUANTITY<float>,float>(std::istream& input,float& d)
{Read_Binary<float>(input,d);}

template<> inline void
Read_Binary<QUANTITY<float>,double>(std::istream& input,double& d)
{Read_Binary<float>(input,d);}

template<> inline void
Read_Binary<QUANTITY<double>,float>(std::istream& input,float& d)
{Read_Binary<double>(input,d);}

template<> inline void
Read_Binary<QUANTITY<double>,double>(std::istream& input,double& d)
{Read_Binary<double>(input,d);}

template<class T> inline void
Read_Primitive(std::ostream& input,QUANTITY<T>& d)
{Read_Primitive(input,d.value);d.unit=UNITS::PARTIAL_UNIT::New();}

template<class T> inline void
Write_Primitive(std::ostream& output,const QUANTITY<T> d)
{Write_Primitive(output,d.value);}

template<class T> inline std::istream& operator>>(std::istream& input,QUANTITY<T>& a)
{input>>a.value;a.unit=UNITS::PARTIAL_UNIT::New();return input;}

template<class T> inline std::ostream& operator<<(std::ostream& output,const QUANTITY<T> a)
{output<<a.value;return output;}

//#####################################################################
}
#endif
