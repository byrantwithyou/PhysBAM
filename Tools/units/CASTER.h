//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Struct CASTER
//#####################################################################
#ifndef __CASTER__
#define __CASTER__

#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include "QUANTITY_FORWARD.h"
namespace PhysBAM{
namespace UNITS{

template<class I>
struct INT_CASTER
{
    I value;

    template<class T>
    explicit INT_CASTER(const T& a)
        :value((I)a)
    {}

    template<class T>
    explicit INT_CASTER(const QUANTITY<T> a)
        :value((I)a.value)
    {
        Unify_One(a.unit);
    }

    operator I() const
    {return value;}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Binary<RW>(output,value);}
};

inline int min(const int a,const INT_CASTER<int> b)
{return PhysBAM::min(a,b.value);}

inline int min(const INT_CASTER<int> a,const int b)
{return PhysBAM::min(a.value,b);}

inline int max(const int a,const INT_CASTER<int> b)
{return PhysBAM::max(a,b.value);}

inline int max(const INT_CASTER<int> a,const int b)
{return PhysBAM::max(a.value,b);}

inline int clamp(const INT_CASTER<int> x,const int xmin,const int xmax)
{return PhysBAM::clamp(x.value,xmin,xmax);}

template<class T> struct CASTER{typedef T TYPE;};
template<> struct CASTER<int>{typedef INT_CASTER<int> TYPE;};
template<> struct CASTER<unsigned char>{typedef INT_CASTER<unsigned char> TYPE;};
template<> struct CASTER<unsigned short>{typedef INT_CASTER<unsigned short> TYPE;};

}
}
#endif
