//#####################################################################
// Copyright 2014, Alexey Stomakhin
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLIP_TRANSLATION_PATH__
#define __FLIP_TRANSLATION_PATH__

#include <Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV>
struct FLIP_TRANSLATION_PATH
{
    typedef typename TV::SCALAR T;

    virtual ~FLIP_TRANSLATION_PATH(){}
    virtual TV X(const T t) const=0;
    virtual TV V(const T t) const=0;
};

template<class TV>
struct FLIP_TRANSLATION_PATH_UNIFORM: public FLIP_TRANSLATION_PATH<TV>
{
    typedef typename TV::SCALAR T;

    TV x;
    TV v;

    FLIP_TRANSLATION_PATH_UNIFORM(const TV& x_input,const TV& v_input):x(x_input),v(v_input){}
    virtual ~FLIP_TRANSLATION_PATH_UNIFORM(){}

    virtual TV X(const T t) const {return x+v*t;}
    virtual TV V(const T t) const {return v;}
};

template<class TV>
struct FLIP_TRANSLATION_PATH_STATIC: public FLIP_TRANSLATION_PATH<TV>
{
    typedef typename TV::SCALAR T;

    TV x;

    FLIP_TRANSLATION_PATH_STATIC(const TV& x_input):x(x_input){}
    virtual ~FLIP_TRANSLATION_PATH_STATIC(){}

    virtual TV X(const T t) const {return x;}
    virtual TV V(const T t) const {return TV();}
};
}
#endif
