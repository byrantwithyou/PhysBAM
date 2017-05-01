//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_OBJECT_EULER  
//#####################################################################
#ifndef __BOUNDARY_OBJECT_EULER__
#define __BOUNDARY_OBJECT_EULER__

#include <Compressible/Conservation_Law_Solvers/BOUNDARY_OBJECT.h>

namespace PhysBAM{

template<class TV> struct GRID_ARRAYS_POLICY;

template<class TV>
class BOUNDARY_OBJECT_EULER:public BOUNDARY_OBJECT<TV,VECTOR<typename TV::SCALAR,TV::m+2> >
{
    typedef typename TV::SCALAR T;typedef VECTOR<T,TV::m+2> TV_DIMENSION;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_DIMENSION_SCALAR::ELEMENT T_ARRAYS_ELEMENT;
    enum{d=TV::m+2};
public:

    BOUNDARY_OBJECT_EULER() = default;
    ~BOUNDARY_OBJECT_EULER() = default;
//#####################################################################
    void Apply_Neumann_Boundary_Condition(TV_DIMENSION& u_1d,const T neumann_face_velocity,const int axis) override;
    void Apply_Neumann_Boundary_Condition(T_ARRAYS_ELEMENT& u_1d,const TV& normal,const T object_velocity_normal_component) override;
//#####################################################################
};
}
#endif
