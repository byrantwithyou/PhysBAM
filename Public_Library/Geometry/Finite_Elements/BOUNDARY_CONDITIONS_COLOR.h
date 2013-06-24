//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CONDITIONS_COLOR
//#####################################################################
#ifndef __BOUNDARY_CONDITIONS_COLOR__
#define __BOUNDARY_CONDITIONS_COLOR__

#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV>
struct BOUNDARY_CONDITIONS_COLOR: public NONCOPYABLE
{
    bool use_discontinuous_velocity;

    BOUNDARY_CONDITIONS_COLOR()
        :use_discontinuous_velocity(false) 
    {}

    virtual ~BOUNDARY_CONDITIONS_COLOR(){}

    virtual TV u_jump(const TV& X,int color0,int color1)=0;
    virtual TV j_surface(const TV& X,int color0,int color1)=0;
};
}
#endif
