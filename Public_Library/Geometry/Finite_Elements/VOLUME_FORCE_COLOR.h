//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VOLUME_FORCE_COLOR
//#####################################################################
#ifndef __VOLUME_FORCE_COLOR__
#define __VOLUME_FORCE_COLOR__

#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV>
struct VOLUME_FORCE_COLOR: public NONCOPYABLE
{
    VOLUME_FORCE_COLOR(){}
    virtual ~VOLUME_FORCE_COLOR(){}

    virtual TV F(const TV& X,int color)=0;
};
}
#endif
