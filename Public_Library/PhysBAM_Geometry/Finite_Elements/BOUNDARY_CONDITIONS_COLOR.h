//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CONDITIONS_COLOR
//#####################################################################
#ifndef __BOUNDARY_CONDITIONS_COLOR__
#define __BOUNDARY_CONDITIONS_COLOR__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV>
struct BOUNDARY_CONDITIONS_COLOR: public NONCOPYABLE
{
    BOUNDARY_CONDITIONS_COLOR(){}
    virtual ~BOUNDARY_CONDITIONS_COLOR(){}

    virtual TV j_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
    virtual TV n_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
    virtual TV d_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
};
}
#endif
