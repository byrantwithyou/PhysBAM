//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTIAL_UNIT
//#####################################################################
#ifndef __PARTIAL_UNIT__
#define __PARTIAL_UNIT__

#include "SIMPLE_UNIT.h"
namespace PhysBAM{
namespace UNITS{

typedef unsigned int VARIABLE;

class PARTIAL_UNIT
{
public:
    SIMPLE_UNIT a; // unit is a x^power
    VARIABLE x;
    rational<int> power;

    PARTIAL_UNIT(const SIMPLE_UNIT& u)
        :a(u),x(0),power(0)
    {}

    PARTIAL_UNIT(const SIMPLE_UNIT& a,const VARIABLE x,const rational<int>& power)
        :a(a),x(x),power(power)
    {}

    static PARTIAL_UNIT Zero()
    {return SIMPLE_UNIT::Zero();}

    static PARTIAL_UNIT One()
    {return SIMPLE_UNIT::One();}

    PARTIAL_UNIT operator*(const PARTIAL_UNIT& u) const
    {Find();u.Find();
    SIMPLE_UNIT ra=a*u.a;
    if(ra.Is_Degenerate()) return ra;
    else if(power==0) return PARTIAL_UNIT(ra,u.x,u.power);
    else if(u.power==0) return PARTIAL_UNIT(ra,x,power);
    else if(x==u.x) return PARTIAL_UNIT(ra,x,power+u.power);
    else return New();}

    PARTIAL_UNIT& operator*=(const PARTIAL_UNIT& u)
    {return *this=*this*u;}

    PARTIAL_UNIT operator/(const PARTIAL_UNIT& u) const
    {return *this*u.Inverse();}

    PARTIAL_UNIT& operator/=(const PARTIAL_UNIT& u)
    {return *this*=u.Inverse();}

    PARTIAL_UNIT Inverse() const
    {return PARTIAL_UNIT(a.Inverse(),x,-power);}

//#####################################################################
    static PARTIAL_UNIT New();
    friend void Unify(const PARTIAL_UNIT& u1,const PARTIAL_UNIT& u2);
    friend void Unify_One(const PARTIAL_UNIT& u);
private:
    void Find() const;
//#####################################################################
};

inline PARTIAL_UNIT operator*(const SIMPLE_UNIT& u1,const PARTIAL_UNIT& u2)
{return PARTIAL_UNIT(u1*u2.a,u2.x,u2.power);}

inline PARTIAL_UNIT
sqrt(const PARTIAL_UNIT& u)
{return PARTIAL_UNIT(sqrt(u.a),u.x,u.power/2);}

inline PARTIAL_UNIT
cbrt(const PARTIAL_UNIT& u)
{return PARTIAL_UNIT(cbrt(u.a),u.x,u.power/3);}

inline PARTIAL_UNIT
pow(const PARTIAL_UNIT& u,const int& y)
{return PARTIAL_UNIT(pow(u.a,y),u.x,u.power*y);}

inline PARTIAL_UNIT
pow(const PARTIAL_UNIT& u,const rational<int>& y)
{return PARTIAL_UNIT(pow(u.a,y),u.x,u.power*y);}

}
}
#endif
