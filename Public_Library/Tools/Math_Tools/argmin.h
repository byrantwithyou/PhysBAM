//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function argmin
//#####################################################################
#ifndef __argmin__
#define __argmin__

namespace PhysBAM{

template<class T>
constexpr inline int argmin(const T a,const T b)
{return a<=b?0:1;}

template<class T>
constexpr inline int argmin(const T a,const T b,const T c)
{return a<=c?(a<=b?0:1):(b<=c?1:2);}

}
#endif
