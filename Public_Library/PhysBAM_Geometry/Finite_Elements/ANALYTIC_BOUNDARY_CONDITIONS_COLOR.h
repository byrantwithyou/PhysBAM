//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_BOUNDARY_CONDITIONS_COLOR
//#####################################################################
#ifndef __ANALYTIC_BOUNDARY_CONDITIONS_COLOR__
#define __ANALYTIC_BOUNDARY_CONDITIONS_COLOR__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_BOUNDARY_CONDITIONS_COLOR: public NONCOPYABLE
{
    typename TV::SCALAR kg,m,s;
    
    ANALYTIC_BOUNDARY_CONDITIONS_COLOR():kg(1),m(1),s(1){}

    virtual TV f_surface(const TV& X,int color1,int color2)=0;
};
}
#endif
