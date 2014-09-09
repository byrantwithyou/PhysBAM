//#####################################################################
// Copyright 2014, Alexey Stomakhin
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_ROTATION_PATH__
#define __MPM_ROTATION_PATH__

#include <Tools/Matrices/ROTATION.h>
#include <Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV>
struct FLIP_ROTATION_PATH
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;

    virtual ~FLIP_ROTATION_PATH(){}
    virtual ROTATION<TV> R(const T t) const=0;
    virtual T_SPIN W(const T t) const=0;
};

template<class TV>
struct FLIP_ROTATION_PATH_UNIFORM:public FLIP_ROTATION_PATH<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;

    const ROTATION<TV> r;
    const T_SPIN w;

    FLIP_ROTATION_PATH_UNIFORM(const ROTATION<TV>& r_input,const T_SPIN& w_input):r(r_input),w(w_input){}
    virtual ~FLIP_ROTATION_PATH_UNIFORM(){}

    virtual ROTATION<TV> R(const T t) const {return ROTATION<TV>::From_Rotation_Vector(w*t)*r;}
    virtual T_SPIN W(const T t) const {return w;}
};

template<class TV>
struct FLIP_ROTATION_PATH_STATIC:public FLIP_ROTATION_PATH<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;

    const ROTATION<TV> r;

    FLIP_ROTATION_PATH_STATIC(const ROTATION<TV>& r_input):r(r_input){}
    virtual ~FLIP_ROTATION_PATH_STATIC(){}

    virtual ROTATION<TV> R(const T t) const {return r;}
    virtual T_SPIN W(const T t) const {return T_SPIN();}
};
}
#endif
