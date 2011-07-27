//#####################################################################
// Copyright 2004, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class JOINT_TYPES_EXAMPLE
//##################################################################### 
#ifndef __JOINT_TYPES_EXAMPLE__
#define __JOINT_TYPES_EXAMPLE__

#include "../ARTICULATED_RIGID_BODIES_3D_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class JOINT_TYPES_EXAMPLE:public ARTICULATED_RIGID_BODIES_3D_EXAMPLE<T,T>
{
public:
    JOINT_TYPES_EXAMPLE(int parameter=0);
    ~JOINT_TYPES_EXAMPLE();
    void Initialize_Rigid_Bodies();
    void Fixed_Joint(const T shift);
    void Ball_Joint(const T shift);
    void Pin_Joint(const T shift);
    void Hinge_Joint(const T shift);
    void Universal_Joint(const T shift);
    void Gimbal_Joint(const T shift);
};
}
#endif
