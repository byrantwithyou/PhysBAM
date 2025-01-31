//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ADVECTION_COLLIDABLE_UNIFORM_FORWARD
//#####################################################################
#ifndef __ADVECTION_COLLIDABLE_UNIFORM_FORWARD__
#define __ADVECTION_COLLIDABLE_UNIFORM_FORWARD__

#include <Incompressible/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_UNIFORM_FORWARD.h>
namespace PhysBAM{

template<class TV,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_UNIFORM<TV> > class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM;
template<class TV,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_BINARY_UNIFORM<TV> > class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_BINARY_UNIFORM;
template<class TV,class T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_UNIFORM<TV> > class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM;
template<class TV,class T_FACE_LOOKUP=FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<TV> > class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM;
}
#endif
