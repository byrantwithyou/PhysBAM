//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Functions Sort and Stable_Sort
//#####################################################################
#ifndef __SORT__
#define __SORT__

#include <algorithm>
namespace PhysBAM{

struct LEXICOGRAPHIC_COMPARE
{
    template<class T_ARRAY> bool operator()(const T_ARRAY& a1,const T_ARRAY& a2) const
    {return std::lexicographical_compare(a1.begin(),a1.end(),a2.begin(),a2.end());}
};
//#####################################################################
}
#endif
