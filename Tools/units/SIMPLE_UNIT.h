//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLE_UNIT
//#####################################################################
#ifndef __SIMPLE_UNIT__
#define __SIMPLE_UNIT__

#include <Tools/Log/DEBUG_UTILITIES.h>
#include <boost/rational.hpp>
namespace PhysBAM{
namespace UNITS{

using boost::rational;

class SIMPLE_UNIT
{
    static const int degenerate=INT_MAX;
public:
    rational<int> m,kg,s; // represents m^m kg^kg s^s

    SIMPLE_UNIT(const rational<int>& m,const rational<int>& kg,const rational<int>& s)
        :m(m),kg(kg),s(s)
    {}

    static SIMPLE_UNIT Zero()
    {return SIMPLE_UNIT(-degenerate,0,0);}

    static SIMPLE_UNIT One()
    {return SIMPLE_UNIT(0,0,0);}

    static SIMPLE_UNIT Infinity()
    {return SIMPLE_UNIT(degenerate,0,0);}

    bool Is_Zero() const
    {return m==-degenerate;}

    bool Is_Infinity() const
    {return m==degenerate;}

    bool Is_Degenerate() const
    {return abs(m)==degenerate;}

    bool operator==(const SIMPLE_UNIT& u) const
    {return m==u.m && kg==u.kg && s==u.s;}

    bool operator!=(const SIMPLE_UNIT& u) const
    {return !(*this==u);}

    SIMPLE_UNIT operator*(const SIMPLE_UNIT& u) const
    {if(Is_Degenerate()){
        if(u.m==-m) PHYSBAM_FATAL_ERROR();
        return *this;}
    if(u.Is_Degenerate()) return u;
    return SIMPLE_UNIT(m+u.m,kg+u.kg,s+u.s);}

    SIMPLE_UNIT operator/(const SIMPLE_UNIT& u) const
    {if(Is_Degenerate()){
        if(u.m==m) PHYSBAM_FATAL_ERROR();
        return *this;}
    if(u.Is_Degenerate()) return u.Inverse();
    return SIMPLE_UNIT(m-u.m,kg-u.kg,s-u.s);}

    SIMPLE_UNIT Inverse() const
    {return SIMPLE_UNIT(-m,-kg,-s);}

//#####################################################################
};

inline SIMPLE_UNIT
sqrt(const SIMPLE_UNIT& u)
{if(u.Is_Degenerate()) return u;
return SIMPLE_UNIT(u.m/2,u.kg/2,u.s/2);}

inline SIMPLE_UNIT
cbrt(const SIMPLE_UNIT& u)
{if(u.Is_Degenerate()) return u;
return SIMPLE_UNIT(u.m/3,u.kg/3,u.s/3);}

inline SIMPLE_UNIT
pow(const SIMPLE_UNIT& u,const int y)
{if(u.Is_Degenerate() && y) return y>0?u:u.Inverse();
return SIMPLE_UNIT(u.m*y,u.kg*y,u.s*y);}

inline SIMPLE_UNIT
pow(const SIMPLE_UNIT& u,const rational<int>& y)
{if(u.Is_Degenerate() && y!=0) return y>0?u:u.Inverse();
return SIMPLE_UNIT(u.m*y,u.kg*y,u.s*y);}

}
}
#endif
