//#####################################################################
// Copyright 2006-2008, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_CONTROL_CALLBACKS
//##################################################################### 
#ifndef __FLUID_CONTROL_CALLBACKS__
#define __FLUID_CONTROL_CALLBACKS__

#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class TV>
class FLUID_CONTROL_CALLBACKS
{
    
    typedef typename TV::SCALAR T;
public:

    FLUID_CONTROL_CALLBACKS()
    {}

    virtual ~FLUID_CONTROL_CALLBACKS()
    {}

//#####################################################################
    virtual T Phi(const TV& location,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return 0;}
    virtual T Face_Velocity(const TV& location,const int axis,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return 0;}
    virtual TV Normal(const TV& location,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return TV();}
    virtual T Shape_Distance_Falloff(const T target_phi,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return 0;}
    virtual T Velocity_Distance_Falloff(const T target_phi,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return 0;}
    virtual T Potential_Distance_Falloff(const T target_phi,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return 0;}
    virtual T Control_Falloff(const TV& location,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return 0;}
//#####################################################################
};
}
#endif
