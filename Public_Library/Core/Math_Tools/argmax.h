//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function argmax
//#####################################################################
#ifndef __argmax__
#define __argmax__

namespace PhysBAM{

template<class T>
constexpr inline int argmax(const T a,const T b)
{return a>=b?0:1;}

template<class T>
constexpr inline int argmax(const T a,const T b,const T c)
{return a>=c?argmax(a,b):(b>=c?1:2);}

}
#endif
