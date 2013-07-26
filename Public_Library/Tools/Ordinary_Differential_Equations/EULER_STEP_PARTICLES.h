//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_STEP_PARTICLES
//#####################################################################
#ifndef __EULER_STEP_PARTICLES__
#define __EULER_STEP_PARTICLES__

#include <Tools/Arrays/ARRAYS_FORWARD.h>

namespace PhysBAM{

template<class TV> struct GRID_ARRAYS_POLICY;

template<class TV>
class EULER_STEP_PARTICLES
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS;
public:
//#####################################################################
    static void Euler_Step_Node(ARRAY_VIEW<TV> X,const GRID<TV>& grid,const ARRAY<TV,TV_INT>& U,const T dt);
    static void Euler_Step_Face(ARRAY_VIEW<TV> X,const GRID<TV>& grid,const T_FACE_ARRAYS& face_velocities,const T dt);
    static void Second_Order_Runge_Kutta_Step(ARRAY_VIEW<TV> X,const GRID<TV>& grid,const ARRAY<TV,TV_INT>& U,const T dt);
//#####################################################################
};
}
#endif
