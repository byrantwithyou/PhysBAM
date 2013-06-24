//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_ON_A_RAY_POLICY 
//#####################################################################
#ifndef __IMPLICIT_OBJECT_ON_A_RAY_POLICY__
#define __IMPLICIT_OBJECT_ON_A_RAY_POLICY__

#include <Tools/Grids_Uniform/UNIFORM_GRID_FORWARD.h>
namespace PhysBAM{

template<class TV> class LEVELSET_IMPLICIT_OBJECT;
template<class T_IMPLICIT_OBJECT> class IMPLICIT_OBJECT_ON_A_RAY;

template<class T_GRID> struct IMPLICIT_OBJECT_ON_A_RAY_POLICY:public IMPLICIT_OBJECT_ON_A_RAY_POLICY<typename T_GRID::GRID_TAG>{};

//#####################################################################
// Uniform
//#####################################################################
template<class TV>
struct IMPLICIT_OBJECT_ON_A_RAY_POLICY<UNIFORM_TAG<TV> >
{
    typedef IMPLICIT_OBJECT_ON_A_RAY<LEVELSET_IMPLICIT_OBJECT<TV> > TYPE;
};
//#####################################################################
}
#endif
