//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INCOMPRESSIBLE_POLICY__
#define __INCOMPRESSIBLE_POLICY__

#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_FORWARD.h>

namespace PhysBAM{

template<class TV> class GRID;

template<class T_GRID> struct INCOMPRESSIBLE_POLICY;

template<class TV> struct INCOMPRESSIBLE_POLICY<GRID<TV> >
{
    typedef INCOMPRESSIBLE_UNIFORM<GRID<TV> > INCOMPRESSIBLE;
    typedef PROJECTION_DYNAMICS_UNIFORM<GRID<TV> > PROJECTION;
};
//#####################################################################
}
#endif
