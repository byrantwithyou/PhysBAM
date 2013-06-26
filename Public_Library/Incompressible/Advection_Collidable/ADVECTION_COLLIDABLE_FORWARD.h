//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ADVECTION_COLLIDABLE_FORWARD
//#####################################################################
#ifndef __ADVECTION_COLLIDABLE_FORWARD__
#define __ADVECTION_COLLIDABLE_FORWARD__

#include <Incompressible/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_UNIFORM_FORWARD.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_NESTED_LOOKUP,class T_NESTED_ADVECTION,class T_FACE_LOOKUP_COLLIDABLE> class ADVECTION_WRAPPER_COLLIDABLE_CELL;
template<class T_GRID,class T2,class T_NESTED_LOOKUP,class T_NESTED_ADVECTION,class T_FACE_LOOKUP_COLLIDABLE> class ADVECTION_WRAPPER_COLLIDABLE_FACE;

}
#endif
