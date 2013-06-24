//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ONE
//#####################################################################
#ifndef __ONE__
#define __ONE__

#include <Tools/Read_Write/READ_WRITE_FUNCTIONS.h>
namespace PhysBAM{

struct ONE
{
    typedef int HAS_UNTYPED_READ_WRITE;

    bool operator!() const
    {return false;}

    ONE operator*(const ONE) const
    {return ONE();}

    bool operator==(const ONE) const
    {return true;}

    ONE Inverse() const
    {return ONE();}

    static ONE One()
    {return ONE();}

    template<class RW> void Read(std::istream& input)
    {}

    template<class RW> void Write(std::ostream& output) const
    {}

//#####################################################################
};

template<class T> inline const T& operator*(const T& x,const ONE)
{return x;}

template<class T> inline const T& operator*(const ONE,const T& x)
{return x;}

template<class T> inline const T& operator/(const T& x,const ONE)
{return x;}

template<class T> inline T& operator*=(T& x,const ONE)
{return x;}

template<class T> inline T& operator/=(T& x,const ONE)
{return x;}

inline std::ostream&
operator<<(std::ostream& output,const ONE)
{return output<<1;}
}
#endif
