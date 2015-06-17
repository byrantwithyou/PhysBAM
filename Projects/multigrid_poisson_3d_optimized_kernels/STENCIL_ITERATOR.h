//#####################################################################
// Copyright 2008-2009, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STENCIL_ITERATOR
//#####################################################################
#ifndef __STENCIL_ITERATOR__
#define __STENCIL_ITERATOR__

#include "STENCIL.h"
namespace PhysBAM{

template<class T,int d>
class STENCIL_ITERATOR
{
    typedef typename conditional<is_const<T>::value,const STENCIL<typename remove_const<T>::type,d>,STENCIL<T,d> >::type T_STENCIL;

    T_STENCIL& stencil;
    int current_index;
public:    

    STENCIL_ITERATOR(T_STENCIL& stencil_input)
        :stencil(stencil_input)
    {Reset();}

    void Reset()
    {current_index=1;}

    const VECTOR<int,d>& Key() const
    {return stencil.entries(current_index).x;}

    T& Data()
    {return stencil.entries(current_index).y;}

    const T& Data() const
    {return stencil.entries(current_index).y;}

    bool Valid() const
    {return current_index<=stencil.Size();}

    void Next()
    {assert(Valid());current_index++;}

//#####################################################################
};
}
#endif
