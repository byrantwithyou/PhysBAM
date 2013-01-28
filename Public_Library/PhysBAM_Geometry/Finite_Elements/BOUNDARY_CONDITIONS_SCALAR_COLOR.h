//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CONDITIONS_SCALAR_COLOR
//#####################################################################
#ifndef __BOUNDARY_CONDITIONS_SCALAR_COLOR__
#define __BOUNDARY_CONDITIONS_SCALAR_COLOR__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV>
struct BOUNDARY_CONDITIONS_SCALAR_COLOR: public NONCOPYABLE
{
    typedef typename TV::SCALAR T;

    bool use_discontinuous_scalar_field;

    BOUNDARY_CONDITIONS_SCALAR_COLOR()
        :use_discontinuous_scalar_field(false) 
    {}

    virtual T u_jump(const TV& X,int color0,int color1)=0;
    virtual T j_surface(const TV& X,int color0,int color1)=0;
};
}
#endif
