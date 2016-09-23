//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Sergey Koltakov, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function max  
//#####################################################################
#ifndef __max__
#define __max__

namespace PhysBAM{

template<class T>
constexpr inline T max(const T a, const T b)
{return (a>b)?a:b;}

template<class T,class ...Args>
constexpr inline T max(const T a,const T b,Args&&... c)
{return max(a,max(b,c...));}

}
#endif
