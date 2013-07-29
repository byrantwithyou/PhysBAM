//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_OBJECT_REFLECTION  
//#####################################################################
#ifndef __BOUNDARY_OBJECT_REFLECTION__
#define __BOUNDARY_OBJECT_REFLECTION__

#include <Compressible/Conservation_Law_Solvers/BOUNDARY_OBJECT.h>

namespace PhysBAM{

template<class TV,class TV_DIMENSION>
class BOUNDARY_OBJECT_REFLECTION:public BOUNDARY_OBJECT<TV,TV_DIMENSION>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_DIMENSION_SCALAR::ELEMENT T_ARRAYS_ELEMENT;
public:
    bool flip;

    BOUNDARY_OBJECT_REFLECTION(const bool flip_input=false)
    {flip=flip_input;}

    ~BOUNDARY_OBJECT_REFLECTION()
    {}

    void Set_Flip(const bool flip_input)
    {flip=flip_input;}

//#####################################################################
    void Apply_Neumann_Boundary_Condition(TV_DIMENSION& u_1d,const T neumann_face_velocity,const int axis) PHYSBAM_OVERRIDE;
    void Apply_Neumann_Boundary_Condition(T_ARRAYS_ELEMENT& u_1d,const TV& normal,const T object_velocity_normal_component) PHYSBAM_OVERRIDE;
//#####################################################################
};
//#####################################################################
// Function Apply_Neumann_Boundary_Condition
//#####################################################################
template<class TV,class TV_DIMENSION> void BOUNDARY_OBJECT_REFLECTION<TV,TV_DIMENSION>::
Apply_Neumann_Boundary_Condition(TV_DIMENSION& u_1d,const T neumann_face_velocity,const int axis)
{
    if(flip) u_1d=u_1d*(T)-1;
}
template<class TV,class TV_DIMENSION> void BOUNDARY_OBJECT_REFLECTION<TV,TV_DIMENSION>::
Apply_Neumann_Boundary_Condition(T_ARRAYS_ELEMENT& u_1d,const TV& normal,const T object_velocity_normal_component)
{
    if(flip) u_1d=u_1d*(T)-1;
}
}
#endif
