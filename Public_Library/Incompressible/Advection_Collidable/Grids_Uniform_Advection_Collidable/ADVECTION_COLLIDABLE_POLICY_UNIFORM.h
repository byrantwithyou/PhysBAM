//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADVECTION_COLLIDABLE_POLICY_UNIFORM__
#define __ADVECTION_COLLIDABLE_POLICY_UNIFORM__

#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_UNIFORM_FORWARD.h>

namespace PhysBAM{

template<class TV>
struct ADVECTION_COLLIDABLE_POLICY
{
private:
    typedef typename TV::SCALAR T;
public:
    // collidable
    typedef ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM<TV,T> ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL;
    typedef ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM<TV> ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE;
    // slip collidable
    typedef ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_SLIP_UNIFORM<TV,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<TV> > ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_SLIP_FACE;
};

}
#endif
