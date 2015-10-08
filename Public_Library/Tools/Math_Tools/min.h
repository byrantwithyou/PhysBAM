//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Sergey Koltakov, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function min 
//#####################################################################
#ifndef __min__
#define __min__

namespace PhysBAM{

template<class T>
constexpr inline T min(const T a, const T b)
{return (a<b)?a:b;}

template<class T,class ...Args>
constexpr inline T min(const T a,const T b,Args&&... c)
{return min(a,min(b,c...));}

}
#endif
