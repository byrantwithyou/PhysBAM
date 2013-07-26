//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_OBJECT_SOLID_VELOCITY  
//#####################################################################
#ifndef __BOUNDARY_OBJECT_SOLID_VELOCITY__
#define __BOUNDARY_OBJECT_SOLID_VELOCITY__

#include <Compressible/Conservation_Law_Solvers/BOUNDARY_OBJECT.h>

namespace PhysBAM{

template<class TV> struct GRID_ARRAYS_POLICY;

template<class TV>
class BOUNDARY_OBJECT_SOLID_VELOCITY:public BOUNDARY_OBJECT<TV,VECTOR<typename TV::SCALAR,TV::m+2> >
{
    typedef typename TV::SCALAR T;typedef VECTOR<T,TV::m+2> TV_DIMENSION;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_DIMENSION_SCALAR::ELEMENT T_ARRAYS_ELEMENT;
    enum{d=TV::m+2};
public:

    BOUNDARY_OBJECT_SOLID_VELOCITY()
    {}

    ~BOUNDARY_OBJECT_SOLID_VELOCITY()
    {}
//#####################################################################
    void Apply_Neumann_Boundary_Condition(TV_DIMENSION& u_1d,const T neumann_face_velocity,const int axis) PHYSBAM_OVERRIDE;
    void Apply_Neumann_Boundary_Condition(T_ARRAYS_ELEMENT& u_1d,const TV& object_velocity,const T unused) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
