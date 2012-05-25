//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_BOUNDARY_CONDITIONS_SCALAR_COLOR
//#####################################################################
#ifndef __ANALYTIC_BOUNDARY_CONDITIONS_SCALAR_COLOR__
#define __ANALYTIC_BOUNDARY_CONDITIONS_SCALAR_COLOR__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_BOUNDARY_CONDITIONS_SCALAR_COLOR: public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    T kg,m,s;
    
    ANALYTIC_BOUNDARY_CONDITIONS_SCALAR_COLOR():kg(1),m(1),s(1){}

    virtual T f_surface(const TV& X,int color0,int color1)=0;
};
}
#endif
