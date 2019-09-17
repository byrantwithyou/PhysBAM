//#####################################################################
// Copyright 2019, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class Make_View
//#####################################################################
#ifndef __Make_View__
#define __Make_View__

#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Matrices/MATRIX_VIEW.h>
namespace PhysBAM{

// Makes a view from a, copying data to tmp if necessary.
// Returns true if a copy was made
template<class T,class T_ARRAY>
enable_if_t<HAS_POINTER<T_ARRAY>::value,bool>
Make_View(ARRAY_VIEW<T>& v,ARRAY_BASE<T,T_ARRAY>& a,ARRAY<T>& tmp)
{
    v.Set(a.Derived());
    return false;
}

template<class T,class T_ARRAY>
enable_if_t<HAS_POINTER<T_ARRAY>::value,bool>
Make_View(ARRAY_VIEW<const T>& v,const ARRAY_BASE<T,T_ARRAY>& a,ARRAY<T>& tmp)
{
    v.Set(a.Derived());
    return false;
}

template<class T,class T_ARRAY>
enable_if_t<!HAS_POINTER<T_ARRAY>::value,bool>
Make_View(ARRAY_VIEW<T>& v,ARRAY_BASE<T,T_ARRAY>& a,ARRAY<T>& tmp)
{
    tmp=a;
    v.Set(tmp);
    return true;
}

template<class T,class T_ARRAY>
enable_if_t<!HAS_POINTER<T_ARRAY>::value,bool>
Make_View(ARRAY_VIEW<const T>& v,const ARRAY_BASE<T,T_ARRAY>& a,ARRAY<T>& tmp)
{
    tmp=a;
    v.Set(tmp);
    return true;
}

}
#endif
