//#####################################################################
// Copyright 2003-2006, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_STATE
//#####################################################################
#ifndef __RIGID_BODY_STATE__
#define __RIGID_BODY_STATE__

#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_MASS.h>
namespace PhysBAM{

template<class TV>
class RIGID_BODY_STATE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    T time;
    FRAME<TV> frame;
    TWIST<TV> twist;
    T_SPIN angular_momentum;

    RIGID_BODY_STATE()
        :time(0)
    {}

    RIGID_BODY_STATE(const FRAME<TV>& frame_input)
        :time(0),frame(frame_input)
    {}

    RIGID_BODY_STATE(const FRAME<TV>& frame_input,const TWIST<TV>& twist_input)
        :time(0),frame(frame_input),twist(twist_input)
    {}

    template<class TV2> explicit RIGID_BODY_STATE(const RIGID_BODY_STATE<TV2>& state_input)
        :time(state_input.time),frame(state_input.frame),twist(state_input.twist),angular_momentum(state_input.angular_momentum)
    {}

    SYMMETRIC_MATRIX<T,0> World_Space_Inertia_Tensor(const DIAGONAL_MATRIX<T,0> moment_of_inertia) const // relative to the center of mass
    {STATIC_ASSERT(TV::m==1);return moment_of_inertia;}

    SYMMETRIC_MATRIX<T,0> World_Space_Inertia_Tensor_Inverse(const DIAGONAL_MATRIX<T,0> moment_of_inertia) const // relative to the center of mass
    {STATIC_ASSERT(TV::m==1);return moment_of_inertia.Inverse();}

    SYMMETRIC_MATRIX<T,1> World_Space_Inertia_Tensor(const DIAGONAL_MATRIX<T,1> moment_of_inertia) const // relative to the center of mass
    {STATIC_ASSERT(TV::m==2);return moment_of_inertia;}

    SYMMETRIC_MATRIX<T,3> World_Space_Inertia_Tensor(const DIAGONAL_MATRIX<T,3>& inertia_tensor) const // relative to the center of mass
    {STATIC_ASSERT(TV::m==3);MATRIX<T,3> object_to_world_transformation=frame.r.Rotation_Matrix();
    return SYMMETRIC_MATRIX<T,3>::Conjugate(object_to_world_transformation,inertia_tensor);}

    RIGID_BODY_MASS<TV,true> World_Space_Rigid_Mass(const RIGID_BODY_MASS<TV>& rigid_mass) const
    {return RIGID_BODY_MASS<TV,true>(rigid_mass.mass,World_Space_Inertia_Tensor(rigid_mass.inertia_tensor));}

    SYMMETRIC_MATRIX<T,1> World_Space_Inertia_Tensor_Inverse(const DIAGONAL_MATRIX<T,1>& moment_of_inertia) const // relative to the center of mass
    {STATIC_ASSERT(TV::m==2);return moment_of_inertia.Inverse();}

    SYMMETRIC_MATRIX<T,3> World_Space_Inertia_Tensor_Inverse(const DIAGONAL_MATRIX<T,3>& inertia_tensor) const // relative to the center of mass
    {STATIC_ASSERT(TV::m==3);MATRIX<T,3> object_to_world_transformation=frame.r.Rotation_Matrix();
    return SYMMETRIC_MATRIX<T,3>::Conjugate(object_to_world_transformation,inertia_tensor.Inverse());}

    RIGID_BODY_MASS<TV,true> World_Space_Rigid_Mass_Inverse(const RIGID_BODY_MASS<TV>& rigid_mass) const
    {return RIGID_BODY_MASS<TV,true>(1/rigid_mass.mass,World_Space_Inertia_Tensor_Inverse(rigid_mass.inertia_tensor));}

    void Update_Angular_Velocity(const DIAGONAL_MATRIX<T,0> inertia_tensor) // needs to be called to keep the angular velocity valid
    {STATIC_ASSERT(TV::m==1);}

    void Update_Angular_Velocity(const DIAGONAL_MATRIX<T,1> moment_of_inertia) // needs to be called to keep the angular velocity valid
    {STATIC_ASSERT(TV::m==2);twist.angular=moment_of_inertia.Solve_Linear_System(angular_momentum);}

    void Update_Angular_Velocity(const DIAGONAL_MATRIX<T,3>& inertia_tensor) // needs to be called to keep the angular velocity valid
    {STATIC_ASSERT(TV::m==3);twist.angular=World_Space_Vector(inertia_tensor.Solve_Linear_System(Object_Space_Vector(angular_momentum)));}

    void Update_Angular_Momentum(const DIAGONAL_MATRIX<T,0> inertia_tensor) // needs to be called to keep the angular velocity valid
    {STATIC_ASSERT(TV::m==1);}

    void Update_Angular_Momentum(const DIAGONAL_MATRIX<T,1> moment_of_inertia) // assumes a valid angular_velocity
    {STATIC_ASSERT(TV::m==2);angular_momentum=moment_of_inertia*twist.angular;}

    void Update_Angular_Momentum(const DIAGONAL_MATRIX<T,3>& inertia_tensor) // assumes a valid angular_velocity
    {STATIC_ASSERT(TV::m==3);angular_momentum=World_Space_Vector(inertia_tensor*Object_Space_Vector(twist.angular));}

    template<class RW> void Read(std::istream& input)
    {char version;Read_Binary<RW>(input,version,time,frame,twist.linear,angular_momentum,twist.angular);assert(version==1);}

    template<class RW> void Write(std::ostream& output) const
    {char version=1;Write_Binary<RW>(output,version,time,frame,twist.linear,angular_momentum,twist.angular);}

    static void Compute_Velocity_Between_States(const RIGID_BODY_STATE& state1,const RIGID_BODY_STATE& state2,RIGID_BODY_STATE& result_state)
    {T one_over_dt=1/(state2.time-state1.time);
    result_state.twist.linear=one_over_dt*(state2.frame.t-state1.frame.t);
    result_state.twist.angular=one_over_dt*(state2.frame.r*state1.frame.r.Inverse()).Rotation_Vector();}

    TV Object_Space_Point(const TV& world_space_point) const
    {return frame.Inverse_Times(world_space_point);}

    TV Object_Space_Vector(const TV& world_space_vector) const
    {return frame.r.Inverse_Rotate(world_space_vector);}

    TV World_Space_Point(const TV& object_space_point) const
    {return frame*object_space_point;}

    TV World_Space_Vector(const TV& object_space_vector) const
    {return frame.r.Rotate(object_space_vector);}

    TV Pointwise_Object_Velocity(const TV& X) const
    {return twist.linear+TV::Cross_Product(twist.angular,X-frame.t);}

//#####################################################################
};
}
#endif
